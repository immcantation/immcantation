#!/usr/bin/env bash
# Script to run a simple pRESTO 0.5.11 workflow on Clontech 5'RACE paired-end data (TakaraBio SMARTer protocol)
#
# Author:  Jason Anthony Vander Heiden
# Date:    2019.05.31
#
# Arguments:
#   -1  Read 1 FASTQ sequence file (sequence beginning with the C-region or J-segment).
#   -2  Read 2 FASTQ sequence file (sequence beginning with the leader or V-segment).
#   -j  C-region reference sequences (reverse complemented).
#       Defaults to /usr/local/share/protocols/Universal/Mouse_TR_CRegion_RC.fasta
#   -r  V-segment reference file.
#       Defaults to /usr/local/share/igblast/fasta/imgt_mouse_tr_v.fasta
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the read 1 filename.
#   -o  Output directory. Will be created if it does not exist.
#       Defaults to a directory matching the sample identifier in the current working directory.
#   -x  The mate-pair coordinate format of the raw data.
#       Defaults to illumina.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -1  Read 1 FASTQ sequence file.\n" \
            "     Sequence beginning with the C-region."
    echo -e "  -2  Read 2 FASTQ sequence file.\n" \
            "     Sequence beginning with the leader."
    echo -e "  -j  C-region reference sequences (reverse complemented).\n" \
            "     Defaults to /usr/local/share/protocols/Universal/Mouse_TR_CRegion_RC.fasta."
    echo -e "  -r  V-segment reference file.\n" \
            "     Defaults to /usr/local/share/igblast/fasta/imgt_mouse_tr_v.fasta."
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the read 1 filename."
    echo -e "  -o  Output directory. Will be created if it does not exist.\n" \
            "     Defaults to a directory matching the sample identifier in the current working directory."
    echo -e "  -x  The mate-pair coordinate format of the raw data.\n" \
            "     Defaults to illumina."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -h  This message."
}

# Argument validation variables
R1_READS_SET=false
R2_READS_SET=false
C_PRIMERS_SET=false
VREF_SEQ_SET=false
OUTNAME_SET=false
OUTDIR_SET=false
NPROC_SET=false
COORD_SET=false

# Get commandline arguments
while getopts "1:2:j:r:n:o:x:p:h" OPT; do
    case "$OPT" in
    1)  R1_READS=$OPTARG
        R1_READS_SET=true
        ;;
    2)  R2_READS=$OPTARG
        R2_READS_SET=true
        ;;
    j)  C_PRIMERS=$OPTARG
        C_PRIMERS_SET=true
        ;;
    r)  VREF_SEQ=$OPTARG
        VREF_SEQ_SET=true
        ;;
    f)  SAMFIELD=$OPTARG
        SAMFIELD_SET=true
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
if ! ${C_PRIMERS_SET}; then
    C_PRIMERS="/usr/local/share/protocols/Universal/Mouse_TR_CRegion_RC.fasta"
elif [ -e ${C_PRIMERS} ]; then
    C_PRIMERS=$(readlink -f ${C_PRIMERS})
else
    echo -e "File ${C_PRIMERS} not found." >&2
    exit 1
fi

# Check reference sequences
if ! ${VREF_SEQ_SET}; then
    VREF_SEQ="/usr/local/share/igblast/fasta/imgt_mouse_tr_v.fasta"
elif [ -e ${VREF_SEQ} ]; then
    VREF_SEQ=$(readlink -f ${VREF_SEQ})
else
    echo -e "File ${VREF_SEQ} not found." >&2
    exit 1
fi

# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true

# AssemblePairs-sequential run parameters
AP_MAXERR=0.3
AP_MINLEN=8
AP_ALPHA=1e-5
AP_MINIDENT=0.5
AP_EVALUE=1e-5
AP_MAXHITS=100

# FilterSeq run parameters
FS_QUAL=20

# MaskPrimers run parameters
MP_MAXERR=0.2
MP_MAXLEN=50
C_FIELD="C_CALL"

# CollapseSeq run parameters
CS_KEEP=true
CS_MISS=0

# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="logs"
PIPELINE_LOG="${LOGDIR}/pipeline-presto.log"
ERROR_LOG="${LOGDIR}/pipeline-presto.err"
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
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
FilterSeq.py quality -s "${OUTNAME}_assemble-pass.fastq" -q $FS_QUAL --nproc $NPROC \
    --outname "${OUTNAME}" --log "${LOGDIR}/quality.log" \
    >> $PIPELINE_LOG  2> $ERROR_LOG
MP_FILE="${OUTNAME}_quality-pass.fastq"
check_error

# Annotate C-region
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s $MP_FILE -p $C_PRIMERS --mode cut --skiprc --revpr\
    --maxlen $MP_MAXLEN --maxerror $MP_MAXERR --pf ${C_FIELD}  \
    --nproc $NPROC --log "${LOGDIR}/cregion.log" --outname "${OUTNAME}" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
CS_FILE="${OUTNAME}_primers-pass.fastq"
check_error

# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
if $CS_KEEP; then
    CollapseSeq.py -s $CS_FILE -n $CS_MISS --inner --keepmiss --uf ${C_FIELD} \
    --outname "${OUTNAME}-final" >> $PIPELINE_LOG 2> $ERROR_LOG
else
    CollapseSeq.py -s $CS_FILE -n $CS_MISS --inner --uf ${C_FIELD} \
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
ParseHeaders.py table -s ${CS_FILE} -f ID ${C_FIELD} \
    --outname "final-total" --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
ParseHeaders.py table -s "${OUTNAME}-final_collapse-unique.fastq" -f ID ${C_FIELD} DUPCOUNT \
    --outname "final-unique" --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
ParseHeaders.py table -s "${OUTNAME}-final_collapse-unique_atleast-2.fastq" -f ID ${C_FIELD} DUPCOUNT \
    --outname "final-unique-atleast2" --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l "${LOGDIR}/assemble.log" \
    -f ID REFID LENGTH OVERLAP GAP ERROR PVALUE EVALUE1 EVALUE2 IDENTITY \
    --outdir ${LOGDIR} > /dev/null 2> $ERROR_LOG &
if $FILTER_LOWQUAL; then
    ParseLog.py -l "${LOGDIR}/quality.log" -f ID QUALITY --outdir ${LOGDIR} \
         > /dev/null 2> $ERROR_LOG &
fi
ParseLog.py -l "${LOGDIR}/cregion.log" -f ID PRSTART PRIMER ERROR \
    --outdir ${LOGDIR} > /dev/null 2> $ERROR_LOG &
wait
check_error

# Zip or delete intermediate and log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compressing files"
LOG_FILES=$(ls ${LOGDIR}/*.log | grep -v "pipeline")
FILTER_FILES="$(basename ${R1_READS})\|$(basename ${R2_READS})\|$(basename ${C_PRIMERS}))"
FILTER_FILES+="\|final_collapse-unique.fastq\|final_collapse-unique_atleast-2.fastq"
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

