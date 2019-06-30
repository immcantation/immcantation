#!/usr/bin/env bash
# Super script to run IgBLAST and Change-O on 10X data
#
# Author:  Jason Anthony Vander Heiden, Ruoyi Jiang
# Date:    2019.05.15
#
# Arguments:
#   -s  FASTA or FASTQ sequence file.
#   -a  10X Cell Ranger V(D)J contig annotation file.
#       Must corresponding with the FASTA/FASTQ input file (all, filtered or consensus).
#   -r  Directory containing IMGT-gapped reference germlines.
#       Defaults to /usr/local/share/germlines/imgt/human/vdj.
#   -g  Species name. One of human or mouse. Defaults to human.
#   -t  Receptor type. One of ig or tr. Defaults to ig.
#   -x  Distance threshold for clonal assignment.
#       If unspecified, clonal assignment is not performed.
#   -m  Distance model for clonal assignment.
#       Defaults to the nucleotide Hamming distance model (ham).
#   -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file
#       and auxillary_data directories. Defaults to /usr/local/share/igblast.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the read 1 filename.
#   -o  Output directory. Will be created if it does not exist.
#       Defaults to a directory matching the sample identifier in the current working directory.
#   -f  Output format. One of changeo or airr. Defaults to changeo.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -i  Specify to allow partial alignments.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -s  FASTA or FASTQ sequence file."
    echo -e "  -a  10X Cell Ranger V(D)J contig annotation CSV file.\n" \
            "     Must corresponding with the FASTA/FASTQ input file (all, filtered or consensus)."
    echo -e "  -r  Directory containing IMGT-gapped reference germlines.\n" \
            "     Defaults to /usr/local/share/germlines/imgt/human/vdj when species is human.\n" \
            "     Defaults to /usr/local/share/germlines/imgt/mouse/vdj when species is mouse."
    echo -e "  -g  Species name. One of human or mouse. Defaults to human."
    echo -e "  -t  Receptor type. One of ig or tr. Defaults to ig."
    echo -e "  -x  Distance threshold for clonal assignment. Specify \"auto\" for automatic detection.\n" \
            "     If unspecified, clonal assignment is not performed."
    echo -e "  -m  Distance model for clonal assignment.\n" \
            "     Defaults to the nucleotide Hamming distance model (ham)."
    echo -e "  -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file\n" \
            "     and auxillary_data directories. Defaults to /usr/local/share/igblast."
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the sequence filename."
    echo -e "  -o  Output directory. Will be created if it does not exist.\n" \
            "     Defaults to a directory matching the sample identifier in the current working directory."
    echo -e "  -f  Output format. One of changeo or airr. Defaults to changeo."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -i  Specify to allow partial alignments."
    echo -e "  -h  This message."
}

# Argument validation variables
READS_SET=false
A10X_SET=false
REFDIR_SET=false
ORGANISM_SET=false
LOCI_SET=false
DIST_SET=false
MODEL_SET=false
IGDATA_SET=false
OUTNAME_SET=false
OUTDIR_SET=false
FORMAT_SET=false
NPROC_SET=false
PARTIAL=""

# Get commandline arguments
while getopts "s:a:r:g:t:x:m:b:n:o:f:p:ih" OPT; do
    case "$OPT" in
    s)  READS=$OPTARG
        READS_SET=true
        ;;
    a)  A10X=$OPTARG
        A10X_SET=true
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
    x)  DIST=$OPTARG
        DIST_SET=true
        ;;
    m)  MODEL=$OPTARG
        MODEL_SET=true
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
    f)  FORMAT=$OPTARG
        FORMAT_SET=true
        ;;
    p)  NPROC=$OPTARG
        NPROC_SET=true
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

if ! ${A10X_SET}; then
    echo -e "You must specify the Cell Ranger annotation file using the -x option." >&2
    exit 1
fi

# Check that files exist and determined absolute paths
if [ -e ${READS} ]; then
    READS=$(readlink -f ${READS})
else
    echo -e "File '${READS}' not found." >&2
    exit 1
fi

if [ -e ${A10X} ]; then
    A10X=$(readlink -f ${A10X})
else
    echo -e "File '${A10X}' not found." >&2
    exit 1
fi

# Set and check species
if ! ${ORGANISM_SET}; then
    ORGANISM="human"
elif [ ${ORGANISM} != "human" ] && [ ${ORGANISM} != "mouse" ]; then
    echo "Species (-g) must be one of 'human' or 'mouse'." >&2
    exit 1
fi

# Set and check receptor type
if ! ${LOCI_SET}; then
    LOCI="ig"
elif [ ${LOCI} != "ig" ] && [ ${LOCI} != "tr" ]; then
    echo "Receptor type (-t) must be one of 'ig' or 'tr'." >&2
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

# Set distance model
if ! ${MODEL_SET}; then
    MODEL="ham"
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

# Check output directory permissions
if [ -e ${OUTDIR} ]; then
    if ! [ -w ${OUTDIR} ]; then
        echo -e "Output directory '${OUTDIR}' is not writable." >&2
        exit 1
    fi
else
    PARENTDIR=$(dirname $(readlink -f ${OUTDIR}))
    if ! [ -w ${PARENTDIR} ]; then
        echo -e "Parent directory '${PARENTDIR}' of new output directory '${OUTDIR}' is not writable." >&2
        exit 1
    fi
fi


# Set format options
if ! ${FORMAT_SET}; then
    FORMAT="changeo"
fi

if [[ "${FORMAT}" == "airr" ]]; then
    EXT="tsv"
    LOCUS_FIELD="locus"
    PROD_FIELD="productive"
else
	EXT="tab"
	LOCUS_FIELD="LOCUS"
	PROD_FIELD="FUNCTIONAL"
fi

# Set number of processes
if ! ${NPROC_SET}; then
    NPROC=$(nproc)
fi

# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true
SPLIT=true
GERMLINES=true
if ! ${DIST_SET}; then
    CLONE=false
else
    CLONE=true
fi

# DefineClones run parameters
DC_MODE="gene"
DC_ACT="set"

# Create germlines parameters
CG_GERM="full dmask"

# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="logs"
PIPELINE_LOG="${LOGDIR}/pipeline-10x.log"
ERROR_LOG="${LOGDIR}/pipeline-10x.err"
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
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Convert to FASTA"
    IG_FILE=$(fastq2fasta.py ${READS})
else
    IG_FILE=${READS}
fi

# Run IgBLAST
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "AssignGenes igblast"
AssignGenes.py igblast -s ${IG_FILE} --organism ${ORGANISM} --loci ${LOCI} \
    -b ${IGDATA} --format blast --nproc ${NPROC} \
    --outname "${OUTNAME}" --outdir . \
     >> $PIPELINE_LOG 2> $ERROR_LOG
FMT7_FILE="${OUTNAME}_igblast.fmt7"
#check_error

# Parse IgBLAST output
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "MakeDb igblast"
MakeDb.py igblast -i ${FMT7_FILE} -s ${IG_FILE} --10x ${A10X} -r ${REFDIR} \
    --extended --failed ${PARTIAL} --outname "${OUTNAME}" --format ${FORMAT} \
    >> $PIPELINE_LOG 2> $ERROR_LOG
DB_PASS="${OUTNAME}_db-pass.${EXT}"
DB_FAIL="${OUTNAME}_db-fail.${EXT}"
check_error

# Split by chain and productivity
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "ParseDb select"
ParseDb.py select -d ${DB_PASS} -f ${LOCUS_FIELD} -u IGH TRB TRD \
    -o "${OUTNAME}_heavy.${EXT}" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
ParseDb.py select -d ${DB_PASS} -f ${LOCUS_FIELD} -u IGK IGL TRA TRG \
    -o "${OUTNAME}_light.${EXT}" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
HEAVY_ALL="${OUTNAME}_heavy.${EXT}"
LIGHT_ALL="${OUTNAME}_light.${EXT}"

printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "ParseDb split"
ParseDb.py split -d "${OUTNAME}_heavy.${EXT}" -f ${PROD_FIELD} \
    >> $PIPELINE_LOG 2> $ERROR_LOG
ParseDb.py split -d "${OUTNAME}_light.${EXT}" -f ${PROD_FIELD} \
    >> $PIPELINE_LOG 2> $ERROR_LOG
HEAVY_PROD="${OUTNAME}_heavy_${PROD_FIELD}-T.${EXT}"
LIGHT_PROD="${OUTNAME}_light_${PROD_FIELD}-T.${EXT}"
HEAVY_NON="${OUTNAME}_heavy_${PROD_FIELD}-F.${EXT}"
LIGHT_NON="${OUTNAME}_light_${PROD_FIELD}-F.${EXT}"

# Assign clones
if $CLONE; then
    if [ "$DIST" == "auto" ]; then
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Detect cloning threshold"
        shazam-threshold -d ${HEAVY_PROD} -m density -n "${OUTNAME}" \
            -f ${FORMAT} -p ${NPROC} \
            > /dev/null 2> $ERROR_LOG
        DIST=$(tail -n1 "${OUTNAME}_threshold-values.tab" | cut -f2)
    fi

    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "DefineClones"
    DefineClones.py -d ${HEAVY_PROD} --model ${MODEL} \
        --dist ${DIST} --mode ${DC_MODE} --act ${DC_ACT} --nproc ${NPROC} \
        --outname "${OUTNAME}_heavy" --log "${LOGDIR}/clone.log" --format ${FORMAT} \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    CLONE_FILE="${OUTNAME}_heavy_clone-pass.${EXT}"
    check_error

    if [ -f "${LIGHT_PROD}" ]; then
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "light_cluster"
        light_cluster.py -d ${CLONE_FILE} -e ${LIGHT_PROD} \
            -o "${OUTNAME}_heavy_clone-light.${EXT}" --format ${FORMAT} --doublets count \
            > /dev/null 2> $ERROR_LOG
        CLONE_FILE="${OUTNAME}_heavy_clone-light.${EXT}"
    else
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Skipping light_cluster"
    fi

    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "CreateGermlines"
    CreateGermlines.py -d ${CLONE_FILE} --cloned -r ${REFDIR} -g ${CG_GERM} \
        --outname "${OUTNAME}_heavy" --log "${LOGDIR}/germline.log" --format ${FORMAT} \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    HEAVY_PROD="${OUTNAME}_heavy_germ-pass.${EXT}"
	check_error
fi

# Zip or delete intermediate files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Compressing files"
TEMP_FILES=$(ls *.tsv *.tab 2>/dev/null | grep -v "${HEAVY_PROD}\|${LIGHT_PROD}\|${HEAVY_NON}\|${LIGHT_NON}")
if [[ ! -z $TEMP_FILES ]]; then
    if $ZIP_FILES; then
        tar -zcf temp_files.tar.gz $TEMP_FILES
    fi
    if $DELETE_FILES; then
        rm $TEMP_FILES
    fi
fi

# End
printf "DONE\n\n"
cd ../
