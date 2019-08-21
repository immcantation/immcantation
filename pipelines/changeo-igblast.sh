#!/usr/bin/env bash
# Super script to run IgBLAST and Change-O
#
# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2018.12.03
#
# Arguments:
#   -s  FASTA or FASTQ sequence file.
#   -r  Directory containing IMGT-gapped reference germlines.
#       Defaults to /usr/local/share/germlines/imgt/human/vdj when species is human.
#       Defaults to /usr/local/share/germlines/imgt/mouse/vdj when species is mouse.
#   -g  Species name. One of human or mouse. Defaults to human.
#   -t  Receptor type. One of ig or tr. Defaults to ig.
#   -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file
#       and auxillary_data directories. Defaults to /usr/local/share/igblast.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the read 1 filename.
#   -o  Output directory. Will be created if it does not exist.
#       Defaults to a directory matching the sample identifier in the current working directory.
#   -f  Output format. One of changeo or airr. Defaults to changeo.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -k  Specify to filter the output to only productive/functional sequences.
#   -i  Specify to allow partial alignments.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -s  FASTA or FASTQ sequence file."
    echo -e "  -r  Directory containing IMGT-gapped reference germlines.\n" \
            "     Defaults to /usr/local/share/germlines/imgt/human/vdj when species is human.\n" \
            "     Defaults to /usr/local/share/germlines/imgt/mouse/vdj when species is mouse."
    echo -e "  -g  Species name. One of human or mouse. Defaults to human."
    echo -e "  -t  Receptor type. One of ig or tr. Defaults to ig."
    echo -e "  -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file\n" \
            "     and auxillary_data directories. Defaults to /usr/local/share/igblast."
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the sequence filename."
    echo -e "  -o  Output directory. Will be created if it does not exist.\n" \
            "     Defaults to a directory matching the sample identifier in the current working directory."
    echo -e "  -f  Output format. One of changeo (default) or airr."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -k  Specify to filter the output to only productive/functional sequences."
    echo -e "  -i  Specify to allow partial alignments."
    echo -e "  -h  This message."
}

# Argument validation variables
READS_SET=false
REFDIR_SET=false
SPECIES_SET=false
RECEPTOR_SET=false
IGDATA_SET=false
OUTNAME_SET=false
OUTDIR_SET=false
FORMAT_SET=false
NPROC_SET=false
FUNCTIONAL=false
PARTIAL=""

# Get commandline arguments
while getopts "s:r:g:t:b:n:o:f:p:kih" OPT; do
    case "$OPT" in
    s)  READS=$OPTARG
        READS_SET=true
        ;;
    r)  REFDIR=$OPTARG
        REFDIR_SET=true
        ;;
    g)  SPECIES=$OPTARG
        SPECIES_SET=true
        ;;
    t)  RECEPTOR=$OPTARG
        RECEPTOR_SET=true
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
    k)  FUNCTIONAL=true
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

# Check that files exist and determined absolute paths
if [ -e ${READS} ]; then
    READS=$(readlink -f ${READS})
else
    echo -e "File '${READS}' not found." >&2
    exit 1
fi

# Set and check species
if ! ${SPECIES_SET}; then
    SPECIES="human"
elif [ ${SPECIES} != "human" ] && [ ${SPECIES} != "mouse" ]; then
    echo "Species (-g) must be one of 'human' or 'mouse'." >&2
    exit 1
fi

# Set and check receptor type
if ! ${RECEPTOR_SET}; then
    RECEPTOR="ig"
elif [ ${RECEPTOR} != "ig" ] && [ ${RECEPTOR} != "tr" ]; then
    echo "Receptor type (-t) must be one of 'ig' or 'tr'." >&2
    exit 1
fi

# Set reference sequence
if ! ${REFDIR_SET}; then
    if [ ${SPECIES} == "human" ]; then
        REFDIR="/usr/local/share/germlines/imgt/human/vdj"
    elif [ ${SPECIES} == "mouse" ]; then
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
GERMLINES=false

# Create germlines parameters
CG_GERM="dmask"

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
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Converting to FASTA"
    IG_FILE=$(fastq2fasta.py ${READS})
else
    IG_FILE=${READS}
fi

# Run IgBLAST
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "IgBLAST"
run_igblast.sh -s ${IG_FILE} -g ${SPECIES} -t ${RECEPTOR} -b ${IGDATA} -n ${NPROC} \
    >> $PIPELINE_LOG 2> $ERROR_LOG
FMT7_FILE=$(basename ${IG_FILE})
FMT7_FILE="${FMT7_FILE%.fasta}.fmt7"
#check_error

# Parse IgBLAST output
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MakeDb igblast"
    MakeDb.py igblast -i ${FMT7_FILE} -s  ${IG_FILE} -r ${REFDIR} \
    --extended --failed ${PARTIAL} \
    --outname "${OUTNAME}" --outdir . --format ${FORMAT} \
    >> $PIPELINE_LOG 2> $ERROR_LOG
    DB_PASS="${OUTNAME}_db-pass.${EXT}"
    DB_FAIL="${OUTNAME}_db-fail.${EXT}"
    LAST_FILE=$DB_PASS
check_error

# Create germlines
if $GERMLINES; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CreateGermlines"
    CreateGermlines.py -d ${LAST_FILE} -r ${REFDIR} -g ${CG_GERM} \
        --outname "${OUTNAME}" --format ${FORMAT} \
        >> $PIPELINE_LOG 2> $ERROR_LOG
	check_error
	GERM_PASS="${OUTNAME}_germ-pass.${EXT}"
	LAST_FILE=$GERM_PASS
fi

if $FUNCTIONAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseDb select"
    ParseDb.py select -d ${LAST_FILE} -f ${PROD_FIELD} -u T TRUE \
        --outname "${OUTNAME}" \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error
    SELECT_PASS="${OUTNAME}_parse-select.${EXT}"
    LAST_FILE=$SELECT_PASS
fi

# Zip or delete intermediate files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compressing files"
TEMP_FILES=$(ls ${DB_PASS} ${DB_FAIL} ${GERM_PASS} ${SELECT_PASS} 2> /dev/null | grep -v "${LAST_FILE}\|$(basename ${READS})")
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
