#!/usr/bin/env bash
# Super script to run IgBLAST and Change-O on 10X data
#
# Author:  Jason Anthony Vander Heiden, Ruoyi Jiang
# Date:    2019.03.14
#
# Arguments:
#   -s  FASTA or FASTQ sequence file.
#   -x  10X Cell Ranger V(D)J contig annotation file.
#       Must corresponding with the FASTA/FASTQ input file (all, filtered or consensus).
#   -r  Directory containing IMGT-gapped reference germlines.
#       Defaults to /usr/local/share/germlines/imgt/human/vdj.
#   -g  Species name. One of human or mouse. Defaults to human.
#   -t  Receptor type. One of ig or tr. Defaults to ig.
#   -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file
#       and auxillary_data directories. Defaults to /usr/local/share/igblast.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the read 1 filename.
#   -o  Output directory.
#       Defaults to the sample name.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -f  Specify to filter the output to only productive/functional sequences.
#   -i  Specify to allow partial alignments.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -s  FASTA or FASTQ sequence file."
    echo -e "  -x  10X Cell Ranger V(D)J contig annotation CSV file.\n" \
            "     Must corresponding with the FASTA/FASTQ input file (all, filtered or consensus)."
    echo -e "  -r  Directory containing IMGT-gapped reference germlines.\n" \
            "     Defaults to /usr/local/share/germlines/imgt/human/vdj when species is human.\n" \
            "     Defaults to /usr/local/share/germlines/imgt/mouse/vdj when species is mouse."
    echo -e "  -g  Species name. One of human or mouse. Defaults to human."
    echo -e "  -t  Receptor type. One of ig or tr. Defaults to ig."
    echo -e "  -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file\n" \
            "     and auxillary_data directories. Defaults to /usr/local/share/igblast."
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the sequence filename."
    echo -e "  -o  Output directory.\n" \
            "     Defaults to the sample name."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -f  Specify to filter the output to only productive/functional sequences."
    echo -e "  -i  Specify to allow partial alignments."
    echo -e "  -h  This message."
}

# Argument validation variables
READS_SET=false
CSV10X_SET=false
REFDIR_SET=false
ORGANISM_SET=false
LOCI_SET=false
IGDATA_SET=false
OUTNAME_SET=false
OUTDIR_SET=false
NPROC_SET=false
FUNCTIONAL=false
PARTIAL=""

# Get commandline arguments
while getopts "s:x:r:g:t:b:n:o:p:fih" OPT; do
    case "$OPT" in
    s)  READS=$OPTARG
        READS_SET=true
        ;;
    x)  CSV10=$OPTARG
        CSV10X_SET=true
        ;;
    r)  REFDIR=$OPTARG
        REFDIR_SET=true
        ;;
    g)  ORGANISM=$OPTARG
        ORGANISM_SET=true
        ;;
    t)  LOCI=$OPTARG
        LOCI_SET=true
        ;;
    b)  IGDATA=$OPTARG
        IGDATA_SET=true
        ;;
    n)  OUTNAME=$OPTARG
        OUTNAME_SET=true
        ;;
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    p)  NPROC=$OPTARG
        NPROC_SET=true
        ;;
    f)  FUNCTIONAL=true
        ;;
    i)  PARTIAL="--partial"
        ;;        
    h)  print_usage
        exit
        ;;
    \?) echo -e "Invalid option: -${OPTARG}" >&2
        exit 1
        ;;
    :)  echo -e "Option -${OPTARG} requires an argument" >&2
        exit 1
        ;;
    esac
done

# Exit if required arguments are not provided
if ! ${READS_SET}; then
    echo -e "You must specify the input sequences using the -s option." >&2
    exit 1
fi

if ! ${CSV10X_SET}; then
    echo -e "You must specify the Cell Ranger annotation file using the -x option." >&2
    exit 1
fi

# Check that files exist and determined absolute paths
if [ -e ${READS} ]; then
    READS=$(readlink -f ${READS})
else
    echo -e "File ${READS} not found." >&2
    exit 1
fi

if [ -e ${CSV10X} ]; then
    CSV10X=$(readlink -f ${CSV10X})
else
    echo -e "File ${CSV10X} not found." >&2
    exit 1
fi

# Set and check species
if ! ${ORGANISM_SET}; then
    ORGANISM="human"
elif [ ${ORGANISM} != "human" ] && [ ${ORGANISM} != "mouse" ]; then
    echo "Species (-g) must be one of human or mouse" >&2
    exit 1
fi

# Set and check receptor type
if ! ${LOCI_SET}; then
    LOCI="ig"
elif [ ${LOCI} != "ig" ] && [ ${LOCI} != "tr" ]; then
    echo "Receptor type (-t) must be one of ig or tr" >&2
    exit 1
fi

# Set reference sequence
if ! ${REFDIR_SET}; then
    if [ ${ORGANISM} == "human" ]; then
        REFDIR="/usr/local/share/germlines/imgt/human/vdj"
    elif [ ${ORGANISM} == "mouse" ]; then
        REFDIR="/usr/local/share/germlines/imgt/mouse/vdj"
    fi
else
    REFDIR=$(readlink -f ${REFDIR})
fi

# Set blast database
if ! ${IGDATA_SET}; then
    IGDATA="/usr/local/share/igblast"
else
    IGDATA=$(readlink -f ${IGDATA})
fi

# Set output name
if ! ${OUTNAME_SET}; then
    OUTNAME=$(basename ${READS} | sed 's/\.[^.]*$//; s/_L[0-9]*_R[0-9]_[0-9]*//')
fi

# Set output directory
if ! ${OUTDIR_SET}; then
    OUTDIR=${OUTNAME}
fi

# Set number of processes
if ! ${NPROC_SET}; then
    NPROC=$(nproc)
fi


# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true
SPLIT=true

# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="logs"
PIPELINE_LOG="${LOGDIR}/pipeline-igblast.log"
ERROR_LOG="${LOGDIR}/pipeline-igblast.err"
mkdir -p ${LOGDIR}
echo '' > $PIPELINE_LOG
echo '' > $ERROR_LOG

# Check for errors
check_error() {
    if [ -s $ERROR_LOG ]; then
        echo -e "ERROR:"
        cat $ERROR_LOG | sed 's/^/    /'
        exit 1
    fi
}

# Set extension
IGBLAST_VERSION=$(igblastn -version  | grep 'Package' |sed s/'Package: '//)
CHANGEO_VERSION=$(python3 -c "import changeo; print('%s-%s' % (changeo.__version__, changeo.__date__))")
if [[ $CHANGEO_VERSION == 0.4* ]]; then
    EXT="tab"
else
	EXT="tab"
fi

# Start
echo -e "IDENTIFIER: ${OUTNAME}"
echo -e "DIRECTORY: ${OUTDIR}"
echo -e "CHANGEO VERSION: ${CHANGEO_VERSION}"
echo -e "IGBLAST VERSION: ${IGBLAST_VERSION}"
echo -e "\nSTART"
STEP=0

# Convert to FASTA if needed
BASE_NAME=$(basename ${READS})
EXT_NAME=${BASE_NAME##*.}
if [ "${EXT_NAME,,}" == "fastq" ] || [ "${EXT_NAME,,}" == "fq" ]; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Convert to FASTA"
    IG_FILE=$(fastq2fasta.py ${READS})
else
    IG_FILE=${READS}
fi

# Run IgBLAST
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssignGenes igblast"
AssignGenes.py igblast -s ${IG_FILE} --organism ${ORGANISM} --loci ${LOCI} -b ${IGDATA}
    --format blast -n ${NPROC} --outname ${OUTNAME} --outdir . \
     >> $PIPELINE_LOG 2> $ERROR_LOG
FMT7_FILE="${OUTNAME}_igblast.fmt7"
#check_error

# Parse IgBLAST output
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MakeDb igblast"
MakeDb.py igblast -i ${FMT7_FILE} -s ${IG_FILE} --10x ${CSV10X} -r ${REFDIR} \
    --scores --regions --failed ${PARTIAL} --outname "${OUTNAME}" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
DB_PASS="${OUTNAME}_db-pass.${EXT}"
DB_FAIL="${OUTNAME}_db-fail.${EXT}"
LAST_FILE=$DB_PASS
check_error

# Split by chain and productivity
if $SPLIT; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseDb split"
    ParseDb.py split -d ${LAST_FILE} -f FUNCTIONAL --outname "${OUTNAME}" \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    ParseDb.py split -d "${OUTNAME}_FUNCTIONAL-T.${EXT}" "${OUTNAME}_FUNCTIONAL-F.${EXT}" -f LOCUS \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    SPLIT_FILES=$(ls "FUNCTIONAL-[TF]_LOCUS-*.${EXT}")
    check_error
fi

# Zip or delete intermediate files
#printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compress files"
#TEMP_FILES=$(ls ${DB_PASS} ${DB_FAIL} 2>/dev/null | grep -v "${LAST_FILE}\|$(basename ${READS})")
#if [[ ! -z $TEMP_FILES ]]; then
#    if $ZIP_FILES; then
#        tar -zcf temp_files.tar.gz $TEMP_FILES
#    fi
#    if $DELETE_FILES; then
#        rm $TEMP_FILES
#    fi
#fi

# End
printf "DONE\n\n"
cd ../
