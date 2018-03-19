#!/usr/bin/env bash
# Super script to run Change-O 0.3.7 cloning and germline reconstruction
#
# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2018.03.19
#
# Arguments:
#   -d  Change-O formatted TSV (TAB) file.
#   -x  Distance threshold for clonal assignment.
#   -r  Directory containing IMGT-gapped reference germlines.
#       Defaults to /usr/local/share/germlines/imgt/human/vdj.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the input filename.
#   -o  Output directory.
#       Defaults to the sample name.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -d  Change-O formatted TSV (TAB) file."
    echo -e "  -x  Distance threshold for clonal assignment."
    echo -e "  -r  Directory containing IMGT-gapped reference germlines.\n" \
            "     Defaults to /usr/local/share/germlines/imgt/human/vdj."
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the input filename."
    echo -e "  -o  Output directory.\n" \
            "     Defaults to the sample name."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -h  This message."
}

# Argument validation variables
DB_SET=false
DIST_SET=false
REFDIR_SET=false
OUTNAME_SET=false
OUTDIR_SET=false
NPROC_SET=false

# Get commandline arguments
while getopts "d:x:r:n:o:p:h" OPT; do
    case "$OPT" in
    d)  DB=${OPTARG}
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
if ! ${DB_SET}; then
    echo -e "You must specify the input database using the -d option." >&2
    exit 1
fi

if ! ${DIST_SET}; then
    echo -e "You must specify the clonal assignment distance threshold using the -x option." >&2
    exit 1
fi

# Set unspecified arguments
if ! ${REFDIR_SET}; then
    REFDIR="/usr/local/share/germlines/imgt/human/vdj"
else
    REFDIR=$(readlink -f ${REFDIR})
fi

if ! ${OUTNAME_SET}; then
    OUTNAME=$(basename ${DB} | sed 's/\.[^.]*$//; s/_L[0-9]*_R[0-9]_[0-9]*//')
fi

if ! ${OUTDIR_SET}; then
    OUTDIR=${OUTNAME}
fi

if ! ${NPROC_SET}; then
    NPROC=$(nproc)
fi

# Check that files exist and determined absolute paths
if [ -e ${DB} ]; then
    DB=$(readlink -f ${DB})
else
    echo -e "File ${DB} not found." >&2
    exit 1
fi

# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true
GERMLINES=true
FUNCTIONAL=true

# DefineClones run parameters
DC_MODEL="ham"
DC_MODE="gene"
DC_ACT="set"

# Create germlines parameters
CG_GERM="dmask"
CG_SFIELD="SEQUENCE_IMGT"
CG_VFIELD="V_CALL"

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

# Set extension
CHANGEO_VERSION=$(python3 -c "import changeo; print('%s-%s' % (changeo.__version__, changeo.__date__))")
if [[ $CHANGEO_VERSION == 0.4* ]]; then
    EXT="tsv"
    DC_COMMAND=""
else
	EXT="tab"
    DC_COMMAND="bygroup"
fi

# Start
echo -e "IDENTIFIER: ${OUTNAME}"
echo -e "DIRECTORY: ${OUTDIR}"
echo -e "CHANGEO VERSION: ${CHANGEO_VERSION}"
echo -e "\nSTART"
STEP=0

if $FUNCTIONAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseDb select"
    ParseDb.py select -d ${DB} -f FUNCTIONAL -u T TRUE \
        --outname "${OUTNAME}" --outdir "${OUTDIR}" \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error
    LAST_FILE="${OUTNAME}_parse-select.${EXT}"
else
    LAST_FILE=${DB}
fi

# Assign clones
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "DefineClones ${DC_COMMAND}"
DefineClones.py ${DC_COMMAND} -d ${LAST_FILE} --model ${DC_MODEL} \
    --dist ${DIST} --mode ${DC_MODE} --act ${DC_ACT} --nproc ${NPROC} \
    --outname "${OUTNAME}" --outdir "${OUTDIR}" --log "${LOGDIR}/clone.log" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
check_error
LAST_FILE="${OUTNAME}_clone-pass.${EXT}"

# Create germlines
if $GERMLINES; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CreateGermlines"
    CreateGermlines.py -d ${LAST_FILE} -r ${REFDIR} -g ${CG_GERM} \
        --sf ${CG_SFIELD} --vf ${CG_VFIELD} --cloned --outname "${OUTNAME}" \
        >> $PIPELINE_LOG 2> $ERROR_LOG
	check_error
	LAST_FILE="${OUTNAME}_germ-pass.${EXT}"
fi

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l "${LOGDIR}/clone.log" -f VALLELE DALLELE JALLELE JUNCLEN SEQUENCES CLONES \
    > /dev/null 2> $ERROR_LOG &
wait

# Zip or delete intermediate and log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compressing files"
LOG_FILES=$(ls ${LOGDIR}/*.log | grep -v "pipeline")
TEMP_FILES=$(ls *.${EXT} | grep -v "${LAST_FILE}\|$(basename ${DB})")
if $ZIP_FILES; then
    tar -zcf log_files.tar.gz $LOG_FILES
    tar -zcf temp_files.tar.gz $TEMP_FILES
fi
if $DELETE_FILES; then
    rm $TEMP_FILES
    rm $LOG_FILES
fi

# End
printf "DONE\n\n"
cd ../