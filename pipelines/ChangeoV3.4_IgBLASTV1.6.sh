#!/usr/bin/env bash
# Super script to run IgBLAST 1.6 and Change-O 0.3.4
#
# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2017.03.30
#
# Arguments:
#   -s  FASTA or FASTQ sequence file.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the read 1 filename.
#   -o  Output directory inside the data directory.
#       Defaults to the sample name.
#   -d  Data directory which serves as the parent of the output directory.
#       Defaults to /data.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -s  FASTA or FASTQ sequence file.\n"
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the sequence filename."
    echo -e "  -o  Output directory inside the data directory.\n" \
            "     Defaults to the sample name."
    echo -e "  -d  Data directory which serves as the parent of the output directory.\n" \
            "     Defaults to /data."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -h  This message."
}

# Argument validation variables
READS=false
OUTNAME_SET=false
OUTDIR_SET=false
DATADIR_SET=false
NPROC_SET=false

# Get commandline arguments
while getopts "s:n:o:d:y:p:h" OPT; do
    case "$OPT" in
    s)  READS=${OPTARG}
        READS_SET=true
        ;;
    n)  OUTNAME=$OPTARG
        OUTNAME_SET=true
        ;;
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    d)  DATADIR=$OPTARG
        DATADIR_SET=true
        ;;
    p)  NPROC=$OPTARG
        NPROC_SET=true
        ;;
    h)  print_usage
        exit
        ;;
    \?) echo -e "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)  echo -e "Option -$OPTARG requires an argument" >&2
        exit 1
        ;;
    esac
done

# Exit if required arguments are not provided
if ! ${READS_SET}; then
    echo -e "You must specify the input sequences using the -s option." >&2
    exit 1
fi

# Set unspecified arguments
if ! ${OUTNAME_SET}; then
    OUTNAME=$(basename ${READS} | sed 's/\.[^.]*$//; s/_L[0-9]*_R[0-9]_[0-9]*//')
fi

if ! ${OUTDIR_SET}; then
    OUTDIR=${OUTNAME}
fi

if ! ${DATADIR_SET}; then
    DATADIR="/data"
fi

if ! ${NPROC_SET}; then
    NPROC=$(nproc)
fi

# Check that files exist and determined absolute paths
if [ -e ${DATADIR}/${READS} ]; then
    READS=$(readlink -f ${DATADIR}/${READS})
else
    echo -e "File ${READS} not found in ${DATADIR}." >&2
    exit 1
fi

# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true
GERMLINES=false
FUNCTIONAL=false

# MakeDb parameters
REF_DIR="/usr/local/share/germlines/imgt/human/vdj"

# Create germlines parameters
CG_GERM="dmask"
CG_SFIELD="SEQUENCE_IMGT"
CG_VFIELD="V_CALL"

# Make output directory
mkdir -p ${DATADIR}/${OUTDIR}; cd ${DATADIR}/${OUTDIR}

# Define log files
LOGDIR="logs"
PIPELINE_LOG="${LOGDIR}/pipeline-igblast.log"
ERROR_LOG="${LOGDIR}/pipeline-igblast.err"
mkdir -p ${LOGDIR}
mkdir -p ${REPORTDIR}
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

# Start
CHANGEO_VERSION=$(python3 -c "import changeo; print('%s-%s' % (changeo.__version__, changeo.__date__))")
IGBLAST_VERSION=$(igblastn -version  | grep 'Package' |sed s/'Package: '//)
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
run_igblast.sh -s ${IG_FILE} -n ${NPROC} \
    >> $PIPELINE_LOG 2> $ERROR_LOG
DB_FILE=$(basename ${IG_FILE})
DB_FILE="${DB_FILE%.fasta}.fmt7"
check_error

# Parse IgBLAST output
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MakeDb igblast"
MakeDb.py igblast -i ${DB_FILE} -s  ${IG_FILE} -r ${REF_DIR} \
    --scores --regions --failed --outname "${OUTNAME}" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
    LAST_FILE="${OUTNAME}_db-pass.tab"
check_error

# Create germlines
if $GERMLINES; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CreateGermlines"
    CreateGermlines.py -d ${LAST_FILE} -r ${REF_DIR} -g ${CG_GERM} \
        --sf ${CG_SFIELD} --vf $CG_VFIELD --outname "${OUTNAME}" \
        >> $PIPELINE_LOG 2> $ERROR_LOG
	check_error
	LAST_FILE="${OUTNAME}_germ-pass.tab"
fi

if $FUNCTIONAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseDb select"
    ParseDb.py select -d ${LAST_FILE} -f FUNCTIONAL -u T --outname "${OUTNAME}" \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error
    LAST_FILE="${OUTNAME}_parse-select.tab"
fi

# Zip or delete intermediate and log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compressing files"
TEMP_FILES=$(ls *.tab | grep -v "${LAST_FILE}")
if $ZIP_FILES; then
    tar -zcf temp_files.tar.gz $TEMP_FILES
fi
if $DELETE_FILES; then
    rm $TEMP_FILES
fi

# End
printf "DONE\n\n"
cd ../