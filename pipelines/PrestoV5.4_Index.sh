#!/usr/bin/env bash
# Script to run the pRESTO 0.5.4 UMI and index correction pipeline on AbVitro AbSeq data
#
# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2017.10.17
#
# Arguments:
#   -x  Identity threshold for clustering of UMI sequences.
#       Defaults to 0.80.
#   -z  Identity threshold for clustering of reads.
#       Defaults to 0.85.
#	-f  Annotation field defining sample name.
#       Defaults to 'SAMPLE'.
#   -o  Output directory.
#       Defaults to 'output'.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS] [FILES...]"
    echo -e "  -x  Identity threshold for clustering of UMI sequences.\n" \
            "     Defaults to 0.8."
    echo -e "  -z  Identity threshold for clustering of reads.\n" \
            "     Defaults to 0.85."
    echo -e "  -f  Annotation field defining sample name.\n" \
            "     Defaults to 'SAMPLE'."
    echo -e "  -o  Output directory.\n" \
            "     Defaults to 'output'."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -h  This message."
}

# Argument validation variables
UMI_IDENT_SET=false
SEQ_IDENT_SET=false
FIELD_SET=false
OUTDIR_SET=false
NPROC_SET=false

# Get commandline arguments
while getopts ":x:z:f:o:p:h" OPT; do
    case "$OPT" in
    x)  UMI_IDENT=$OPTARG
        UMI_IDENT_SET=true
        ;;
    z)  SEQ_IDENT=$OPTARG
        SEQ_IDENT_SET=true
        ;;        
    f)  FIELD=$OPTARG
        FIELD_SET=true
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
    \?) echo -e "Invalid option: -${OPTARG}" >&2
        exit 1
        ;;
    :)  echo -e "Option -${OPTARG} requires an argument" >&2
        exit 1
        ;;
    esac
done
shift "$((OPTIND - 1))"

# Exit if required arguments are not provided
if [ $# -lt 2 ]; then
    echo -e "You must specify a list of sequence files as positional arguments." >&2
    exit 1
else
    READS=$(echo $@ | xargs readlink -f)
fi

# Check reads
for R in $READS; do
	if [ ! -e $R ]; then
		echo -e "File ${R} not found." >&2
		exit 1
	fi	
done

# Set unspecified arguments
if ! ${UMI_IDENT_SET}; then
    UMI_IDENT=0.80
fi

if ! ${SEQ_IDENT_SET}; then
    SEQ_IDENT=0.85
fi

if ! ${FIELD_SET}; then
    FIELD="SAMPLE"
fi

if ! ${OUTDIR_SET}; then
    OUTDIR="output"
fi

if ! ${NPROC_SET}; then
    NPROC=$(nproc)
fi


# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true

# ClusterSets-barcode run parameters
BARCODE_FIELD="BARCODE"

# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="logs"
PIPELINE_LOG="${LOGDIR}/pipeline-indexing.log"
ERROR_LOG="${LOGDIR}/pipeline-indexing.err"
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

# Start
PRESTO_VERSION=$(python3 -c "import presto; print('%s-%s' % (presto.__version__, presto.__date__))")
echo -e "DIRECTORY: ${OUTDIR}"
echo -e "PRESTO VERSION: ${PRESTO_VERSION}"
echo -e "\nSTART"
STEP=0

# Cluster UMIs
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Merging files"
merge_fastq.py merged.fastq $READS

# Cluster UMIs
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ClusterSets barcode"
ClusterSets.py barcode -s merged.fastq -f $BARCODE_FIELD -k CLUSTER --ident $UMI_IDENT --prefix UMI \
    --nproc $NPROC --log "${LOGDIR}/cluster-barcode.log" --outdir . \
    >> $PIPELINE_LOG 2> $ERROR_LOG
#check_error

# Cluster reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ClusterSets set"
ClusterSets.py set -s merged_cluster-pass.fastq -f CLUSTER -k CLUSTER --ident $SEQ_IDENT --prefix SEQ \
    --nproc $NPROC --log "${LOGDIR}/cluster-set.log" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
#check_error

# Generate consensus cluster identifiers
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders collapse"
ParseHeaders.py collapse -s merged_cluster-pass_cluster-pass.fastq -f CLUSTER --act cat \
	>> $PIPELINE_LOG 2> $ERROR_LOG
#check_error
	
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "UnifyHeaders consensus"
UnifyHeaders.py consensus -s merged_cluster-pass_cluster-pass_reheader.fastq \
	-f CLUSTER -k $FIELD --log "${LOGDIR}/unify.log" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
#check_error

# Split files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s merged_cluster-pass_cluster-pass_reheader_unify-pass.fastq \
	-f $FIELD --outname final >> $PIPELINE_LOG 2> $ERROR_LOG
#check_error

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l "${LOGDIR}/cluster-barcode.log" "${LOGDIR}/cluster-set.log" \
    -f ID \
    --outdir ${LOGDIR} > /dev/null 2> $ERROR_LOG &
ParseLog.py -l "${LOGDIR}/unify.log" -f ID \
    --outdir ${LOGDIR} > /dev/null 2> $ERROR_LOG &
wait
#check_error

# Zip or delete intermediate and log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compressing files"
LOG_FILES=$(ls ${LOGDIR}/*.log | grep -v "pipeline")
FILTER_FILES="final_"
TEMP_FILES=$(ls *.fastq | grep -v ${FILTER_FILES})
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

