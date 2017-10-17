#!/usr/bin/env bash
# Script to run the pRESTO 0.5.4 initial assembly and annotation on AbVitro AbSeq data
#
# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2017.10.16
#
# Arguments:
#   -1  Read 1 FASTQ sequence file (sequence beginning with the C-region or J-segment).
#   -2  Read 2 FASTQ sequence file (sequence beginning with the leader or V-segment).
#   -j  Read 1 FASTA primer sequences (C-region or J-segment).
#       Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R1_Human_IG_Primers.fasta
#   -v  Read 2 FASTA primer sequences (template switch or V-segment).
#       Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R2_TS.fasta.
#   -r  V-segment reference file.
#       Defaults to /usr/local/share/germlines/igblast/fasta/imgt_human_ig_v.fasta
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the read 1 filename.
#   -o  Output directory.
#       Defaults to the sample name.
#   -x  The mate-pair coordinate format of the raw data.
#       Defaults to illumina.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -1  Read 1 FASTQ sequence file.\n" \
            "     Sequence beginning with the C-region or J-segment)."
    echo -e "  -2  Read 2 FASTQ sequence file.\n" \
            "     Sequence beginning with the leader or V-segment)."
    echo -e "  -j  Read 1 FASTA primer sequences.\n" \
            "     Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R1_Human_IG_Primers.fasta."
    echo -e "  -v  Read 2 FASTA primer or template switch sequences.\n" \
            "     Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R2_TS.fasta."
    echo -e "  -r  V-segment reference file.\n" \
            "     Defaults to /usr/local/share/igblast/fasta/imgt_human_ig_v.fasta."
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the read 1 filename."
    echo -e "  -o  Output directory.\n" \
            "     Defaults to the sample name."
    echo -e "  -x  The mate-pair coordinate format of the raw data.\n" \
            "     Defaults to illumina."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -h  This message."
}

# Argument validation variables
R1_READS_SET=false
R2_READS_SET=false
R1_PRIMERS_SET=false
R2_PRIMERS_SET=false
VREF_SEQ_SET=false
OUTNAME_SET=false
OUTDIR_SET=false
NPROC_SET=false
COORD_SET=false

# Get commandline arguments
while getopts "1:2:j:v:r:n:o:x:p:h" OPT; do
    case "$OPT" in
    1)  R1_READS=${OPTARG}
        R1_READS_SET=true
        ;;
    2)  R2_READS=${OPTARG}
        R2_READS_SET=true
        ;;
    j)  R1_PRIMERS=${OPTARG}
        R1_PRIMERS_SET=true
        ;;
    v)  R2_PRIMERS=${OPTARG}
        R2_PRIMERS_SET=true
        ;;
    r)  VREF_SEQ=${OPTARG}
        VREF_SEQ_SET=true
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
    \?) echo -e "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)  echo -e "Option -$OPTARG requires an argument" >&2
        exit 1
        ;;
    esac
done

# Exit if required arguments are not provided
if ! ${R1_READS_SET} || ! ${R2_READS_SET}; then
    echo -e "You must specify both read files using the -1 and -2 options." >&2
    exit 1
fi

# Set unspecified arguments
if ! ${OUTNAME_SET}; then
    OUTNAME=$(basename ${R1_READS} | sed 's/\.[^.]*$//; s/_L[0-9]*_R[0-9]_[0-9]*//')
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
if [ -e ${R1_READS} ]; then
    R1_READS=$(readlink -f ${R1_READS})
else
    echo -e "File ${R1_READS} not found." >&2
    exit 1
fi

# Check R2 reads
if [ -e ${R2_READS} ]; then
    R2_READS=$(readlink -f ${R2_READS})
else
    echo -e "File ${R2_READS} not found." >&2
    exit 1
fi

# Check R1 primers
if ! ${R1_PRIMERS_SET}; then
    R1_PRIMERS="/usr/local/share/protocols/AbSeq/AbSeq_R1_Human_IG_Primers.fasta"
elif [ -e ${R1_PRIMERS} ]; then
    R1_PRIMERS=$(readlink -f ${R1_PRIMERS})
else
    echo -e "File ${R1_PRIMERS} not found." >&2
    exit 1
fi

# Check R2 primers
if ! ${R2_PRIMERS_SET}; then
    R2_PRIMERS="/usr/local/share/protocols/AbSeq/AbSeq_R2_TS.fasta"
elif [ -e ${R2_PRIMERS} ]; then
    R2_PRIMERS=$(readlink -f ${R2_PRIMERS})
else
    echo -e "File ${R2_PRIMERS} not found." >&2
    exit 1
fi

# Check reference sequences
if ! ${VREF_SEQ_SET}; then
    VREF_SEQ="/usr/local/share/igblast/fasta/imgt_human_ig_v.fasta"
elif [ -e ${VREF_SEQ} ]; then
    VREF_SEQ=$(readlink -f ${VREF_SEQ})
else
    echo -e "File ${VREF_SEQ} not found." >&2
    exit 1
fi

# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true
FILTER_LOWQUAL=true

# AssemblePairs-sequential run parameters
AP_MAXERR=0.3
AP_MINLEN=8
AP_ALPHA=1e-5
AP_MINIDENT=0.5
AP_EVALUE=1e-5
AP_MAXHITS=100

# FilterSeq run parameters
FS_QUAL=20
FS_MASK=30

# MaskPrimers run parameters
MP_UIDLEN=17
MP_R1_MAXERR=0.2
MP_R2_MAXERR=0.5

# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="logs"
PIPELINE_LOG="${LOGDIR}/pipeline-assemble.log"
ERROR_LOG="${LOGDIR}/pipeline-assemble.err"
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
echo -e "IDENTIFIER: ${OUTNAME}"
echo -e "DIRECTORY: ${OUTDIR}"
echo -e "PRESTO VERSION: ${PRESTO_VERSION}"
echo -e "\nSTART"
STEP=0

# Assemble paired ends via mate-pair alignment
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs sequential"
AssemblePairs.py sequential -1 $R2_READS -2 $R1_READS -r $VREF_SEQ \
    --coord $COORD --rc tail --minlen $AP_MINLEN --maxerror $AP_MAXERR --alpha $AP_ALPHA \
    --scanrev --minident $AP_MINIDENT --evalue $AP_EVALUE --maxhits $AP_MAXHITS --aligner blastn \
    --nproc $NPROC --log "${LOGDIR}/assemble.log" --outname "${OUTNAME}" --outdir . \
    >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Remove low quality reads
if $FILTER_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
    #OUTPREFIX="$(printf '%02d' $STEP)--${OUTNAME}"
    FilterSeq.py quality -s "${OUTNAME}_assemble-pass.fastq" -q $FS_QUAL --nproc $NPROC \
        --outname "${OUTNAME}" --log "${LOGDIR}/quality.log" \
        >> $PIPELINE_LOG  2> $ERROR_LOG
    MP_FILE="${OUTNAME}_quality-pass.fastq"
    check_error
else
    MP_FILE="${OUTNAME}_assemble-pass.fastq"
fi

# Identify primers and UID
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
MaskPrimers.py score -s $MP_FILE -p $R2_PRIMERS --mode cut \
    --start $MP_UIDLEN --barcode --maxerror $MP_R2_MAXERR --nproc $NPROC \
    --log "${LOGDIR}/primers-2.log" --outname "${OUTNAME}" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
MaskPrimers.py score -s "${OUTNAME}_primers-pass.fastq" -p $R1_PRIMERS --mode cut \
    --start 0 --maxerror $MP_R1_MAXERR --revpr --nproc $NPROC \
    --log "${LOGDIR}/primers-1.log" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Add sample annotation
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders"
ParseHeaders.py add -s "${OUTNAME}_primers-pass_primers-pass.fastq" -f SAMPLE -u $OUTNAME \
    --outname "${OUTNAME}"  > /dev/null 2> $ERROR_LOG &
check_error

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l "${LOGDIR}/assemble.log" \
    -f ID REFID LENGTH OVERLAP GAP ERROR PVALUE EVALUE1 EVALUE2 IDENTITY FIELDS1 FIELDS2 \
    --outdir ${LOGDIR} > /dev/null 2> $ERROR_LOG &
if $FILTER_LOWQUAL; then
    ParseLog.py -l "${LOGDIR}/quality.log" -f ID QUALITY --outdir ${LOGDIR} \
         > /dev/null 2> $ERROR_LOG &
fi
ParseLog.py -l "${LOGDIR}/primers-1.log" "${LOGDIR}/primers-2.log" -f ID BARCODE PRIMER ERROR \
    --outdir ${LOGDIR} > /dev/null 2> $ERROR_LOG &
wait
check_error

# Zip or delete intermediate and log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compressing files"
LOG_FILES=$(ls ${LOGDIR}/*.log | grep -v "pipeline")
FILTER_FILES="$(basename R1_READS)\|$(basename R2_READS)\|$(basename R1_PRIMERS)\|$(basename R2_PRIMERS)"
FILTER_FILES+="\|reheader.fastq"
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

