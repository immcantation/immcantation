#!/usr/bin/env bash
# Super script to run Change-O 0.3.7 cloning and germline reconstruction
#
# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2018.07.03
#
# Arguments:
#   -d  Change-O formatted TSV (TAB) file.
#   -x  Distance threshold for clonal assignment.
#   -m  Distance model for clonal assignment.
#       Defaults to the nucleotide Hamming distance model (ham).
#   -r  Directory containing IMGT-gapped reference germlines.
#       Defaults to /usr/local/share/germlines/imgt/human/vdj.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the input filename.
#   -o  Output directory. Will be created if it does not exist.
#       Defaults to a directory matching the sample identifier in the current working directory.
#   -f  Output format. One of changeo or airr. Defaults to changeo.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -a  Specify to clone the full data set.
#       By default the data will be filtering to only productive/functional sequences.
#   -z  Specify to disable cleaning and compression of temporary files.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -d  Change-O formatted TSV (TAB) file."
    echo -e "  -x  Distance threshold for clonal assignment."
    echo -e "  -m  Distance model for clonal assignment.\n" \
            "     Defaults to the nucleotide Hamming distance model (ham)."
    echo -e "  -r  Directory containing IMGT-gapped reference germlines.\n" \
            "     Defaults to /usr/local/share/germlines/imgt/human/vdj."
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the input filename."
    echo -e "  -o  Output directory. Will be created if it does not exist.\n" \
            "     Defaults to a directory matching the sample identifier in the current working directory."
    echo -e "  -f  Output format. One of airr (default) or changeo."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -a  Specify to clone the full data set.\n" \
            "     By default the data will be filtering to only productive/functional sequences."
    echo -e "  -z  Specify to disable cleaning and compression of temporary files."
    echo -e "  -h  This message."
}

# Argument validation variables
DB_SET=false
DIST_SET=false
MODEL_SET=false
REFDIR_SET=false
OUTNAME_SET=false
OUTDIR_SET=false
FORMAT_SET=false
NPROC_SET=false

# Argument defaults
FUNCTIONAL=true
ZIP_FILES=true
DELETE_FILES=true

# Get commandline arguments
while getopts "d:x:m:r:n:o:f:p:azh" OPT; do
    case "$OPT" in
    d)  DB=$OPTARG
        DB_SET=true
        ;;
    x)  DIST=$OPTARG
        DIST_SET=true
        ;;
    r)  REFDIR=$OPTARG
        REFDIR_SET=true
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
    a)  FUNCTIONAL=false
        ;;
    z)  ZIP_FILES=false
        DELETE_FILES=false
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
if ! ${DB_SET}; then
    echo -e "You must specify the input database using the -d option." >&2
    exit 1
fi

if ! ${DIST_SET}; then
    echo -e "You must specify the clonal assignment distance threshold using the -x option." >&2
    exit 1
fi

# Set unspecified arguments
if ! ${MODEL_SET}; then
    MODEL="ham"
fi

if ! ${REFDIR_SET}; then
    REFDIR="/usr/local/share/germlines/imgt/human/vdj"
else
    REFDIR=$(realpath ${REFDIR})
fi

if ! ${OUTNAME_SET}; then
    OUTNAME=$(basename ${DB} | sed 's/\.[^.]*$//; s/_L[0-9]*_R[0-9]_[0-9]*//')
fi

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
    PARENTDIR=$(dirname $(realpath ${OUTDIR}))
    if ! [ -w ${PARENTDIR} ]; then
        echo -e "Parent directory '${PARENTDIR}' of new output directory '${OUTDIR}' is not writable." >&2
        exit 1
    fi
fi

# Set format options
if ! ${FORMAT_SET}; then
    FORMAT="airr"
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

# Process settings
if ! ${NPROC_SET}; then
    NPROC=$(nproc)
fi

# Check that files exist and determined absolute paths
if [ -e ${DB} ]; then
    DB=$(realpath ${DB})
else
    echo -e "File ${DB} not found." >&2
    exit 1
fi

# Define pipeline steps
GERMLINES=true

# DefineClones run parameters
DC_MODE="gene"
DC_ACT="set"

# Create germlines parameters
CG_GERM="dmask"

# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="logs"
PIPELINE_LOG="${LOGDIR}/pipeline-clone.log"
ERROR_LOG="${LOGDIR}/pipeline-clone.err"
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

# Get changeo version
CHANGEO_VERSION=$(python3 -c "import changeo; print('%s-%s' % (changeo.__version__, changeo.__date__))")

# Start
echo -e "IDENTIFIER: ${OUTNAME}"
echo -e "DIRECTORY: ${OUTDIR}"
echo -e "CHANGEO VERSION: ${CHANGEO_VERSION}"
echo -e "\nSTART"
STEP=0

if $FUNCTIONAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseDb select"
    ParseDb.py select -d ${DB} -f ${PROD_FIELD} -u T TRUE \
        --outname "${OUTNAME}" --outdir . \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error
    SELECT_PASS="${OUTNAME}_parse-select.${EXT}"
    LAST_FILE=$SELECT_PASS
else
    LAST_FILE=${DB}
fi

# Assign clones
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "DefineClones"
DefineClones.py -d ${LAST_FILE} --model ${MODEL} \
    --dist ${DIST} --mode ${DC_MODE} --act ${DC_ACT} --nproc ${NPROC} \
    --outname "${OUTNAME}" --outdir . --format ${FORMAT} --log "${LOGDIR}/clone.log" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
check_error
CLONE_PASS="${OUTNAME}_clone-pass.${EXT}"
LAST_FILE=$CLONE_PASS

# Create germlines
if $GERMLINES; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CreateGermlines"
    CreateGermlines.py -d ${LAST_FILE} -r ${REFDIR} -g ${CG_GERM} \
        --cloned --outname "${OUTNAME}" --format ${FORMAT} \
        >> $PIPELINE_LOG 2> $ERROR_LOG
	check_error
	GERM_PASS="${OUTNAME}_germ-pass.${EXT}"
	LAST_FILE=$GERM_PASS
fi

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l "${LOGDIR}/clone.log" -f VALLELE DALLELE JALLELE JUNCLEN SEQUENCES CLONES \
    > /dev/null 2> $ERROR_LOG &
wait

# Zip or delete log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compressing files"
LOG_FILES=$(ls ${LOGDIR}/*.log | grep -v "pipeline")
if [[ ! -z $LOG_FILES ]]; then
    if $ZIP_FILES; then
        tar -zcf log_files.tar.gz $LOG_FILES
    fi
    if $DELETE_FILES; then
        rm $LOG_FILES
    fi
fi
# Zip or delete intermediate files
TEMP_FILES=$(ls ${SELECT_PASS} ${CLONE_PASS} ${GERM_PASS}  2> /dev/null | grep -v "${LAST_FILE}\|$(basename ${DB})")
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
