#!/usr/bin/env bash
# Super script to run simple pRESTO filtering and duplicate removal
# 
# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2018.09.30
# 
# Arguments:
#   -s  FASTQ sequence file.
#   -y  YAML file providing description fields for report generation.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the read 1 filename.
#   -o  Output directory. Will be created if it does not exist.
#       Defaults to a directory matching the sample identifier in the current working directory.
#   -x  The header format of the raw data.
#       Defaults to illumina.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -s  FASTQ sequence file."
    echo -e "  -y  YAML file providing description fields for report generation."
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the read 1 filename."
    echo -e "  -o  Output directory. Will be created if it does not exist.\n" \
            "     Defaults to a directory matching the sample identifier in the current working directory."
    echo -e "  -x  The header format of the raw data.\n" \
            "     Defaults to illumina."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -h  This message."
}

# Argument validation variables
READS_SET=false
YAML_SET=FALSE
OUTNAME_SET=false
OUTDIR_SET=false
NPROC_SET=false
COORD_SET=false

# Get commandline arguments
while getopts "s:y:n:o:x:p:h" OPT; do
    case "$OPT" in
    s)  READS=$OPTARG
        READS_SET=true
        ;;
    y)  YAML=$OPTARG
        YAML_SET=true
        ;;
    n)  OUTNAME=$OPTARG
        OUTNAME_SET=true
        ;;
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    x)  COORD=$OPTARG
        COORD_SET=true
        ;;
    p)  NPROC=$OPTARG
        NPROC_SET=true
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
    echo -e "You must specify the read file using the -s option." >&2
    exit 1
fi

if ! ${YAML_SET}; then
    echo -e "You must specify the description file in YAML format using the -y option." >&2
    exit 1
fi

# Set unspecified arguments
if ! ${OUTNAME_SET}; then
    OUTNAME=$(basename ${READS} | sed 's/\.[^.]*$//; s/_L[0-9]*_R[0-9]_[0-9]*//')
fi

if ! ${OUTDIR_SET}; then
    OUTDIR=${OUTNAME}
fi

if ! ${NPROC_SET}; then
    NPROC=$(nproc)
fi

if ! ${COORD_SET}; then
    COORD="illumina"
fi

# Check R1 reads
if [ -e ${READS} ]; then
    READS=$(readlink -f ${READS})
else
    echo -e "File ${READS} not found." >&2
    exit 1
fi

# Check report yaml file
if [ -e ${YAML} ]; then
    YAML=$(readlink -f ${YAML})
else
    echo -e "File ${YAML} not found." >&2
    exit 1
fi

# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true
CONVERT_HEADERS=true
FILTER_LOWQUAL=true
REPORT=true

# FilterSeq run parameters
FS_QUAL=20

# CollapseSeq run parameters
CS_KEEP=false
CS_MISS=20

# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="logs"
REPORTDIR="report"
PIPELINE_LOG="${LOGDIR}/pipeline.log"
ERROR_LOG="${LOGDIR}/pipeline.err"
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
PRESTO_VERSION=$(python3 -c "import presto; print('%s-%s' % (presto.__version__, presto.__date__))")
echo -e "IDENTIFIER: ${OUTNAME}"
echo -e "DIRECTORY: ${OUTDIR}"
echo -e "PRESTO VERSION: ${PRESTO_VERSION}"
echo -e "\nSTART"
STEP=0

# Remove low quality reads
if $CONVERT_HEADERS; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ConvertHeaders ${COORD}"
    ConvertHeaders.py $COORD -s $READS --outname "${OUTNAME}" --outdir . \
        >> $PIPELINE_LOG  2> $ERROR_LOG
    FS_FILE="${OUTNAME}_convert-pass.fastq"
    check_error
else
    FS_FILE=$READS
fi


# Remove low quality reads
if $FILTER_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
    FilterSeq.py quality -s $FS_FILE -q $FS_QUAL --nproc $NPROC \
        --outname "${OUTNAME}" --outdir . --log "${LOGDIR}/quality.log" \
        >> $PIPELINE_LOG  2> $ERROR_LOG
    CS_FILE="${OUTNAME}_quality-pass.fastq"
    check_error
else
    CS_FILE=$FS_FILE
fi


# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
if $CS_KEEP; then
    CollapseSeq.py -s $CS_FILE -n $CS_MISS --inner --keepmiss \
    --outname "${OUTNAME}-final" >> $PIPELINE_LOG 2> $ERROR_LOG
else
    CollapseSeq.py -s $CS_FILE -n $CS_MISS --inner \
    --outname "${OUTNAME}-final" >> $PIPELINE_LOG 2> $ERROR_LOG
fi
check_error


# Filter to sequences with at least 2 supporting sources
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s "${OUTNAME}-final_collapse-unique.fastq" -f DUPCOUNT --num 2 \
    >> $PIPELINE_LOG 2> $ERROR_LOG
check_error


# Create table of final repertoire
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
ParseHeaders.py table -s "${OUTNAME}-final_collapse-unique.fastq" \
    -f ID DUPCOUNT --outname "final-unique" \
    --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
ParseHeaders.py table -s "${OUTNAME}-final_collapse-unique_atleast-2.fastq" \
    -f ID DUPCOUNT --outname "final-unique-atleast2" \
    --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
check_error


# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
if $FILTER_LOWQUAL; then
    ParseLog.py -l "${LOGDIR}/quality.log" -f ID QUALITY \
        --outdir ${LOGDIR} > /dev/null &
fi
wait
check_error


# Generate pRESTO report
if $REPORT; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Generating report"
    REPORT_SCRIPT="report_abseq3(\"${LOGDIR}\", sample=\"${OUTNAME}\", output_dir=\"${REPORTDIR}\", config=\"${YAML}\", quiet=FALSE)"
    Rscript -e "library(prestor); ${REPORT_SCRIPT}" > ${REPORTDIR}/report.out 2> ${REPORTDIR}/report.err
fi


# Zip or delete intermediate and log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compressing files"
LOG_FILES=$(ls ${LOGDIR}/*.log | grep -v "pipeline")
FILTER_FILES="$(basename ${READS})"
FILTER_FILES+="\|final_total.fastq\|final_collapse-unique.fastq\|final_collapse-unique_atleast-2.fastq"
TEMP_FILES=$(ls *.fastq  2>/dev/null | grep -v ${FILTER_FILES})
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

