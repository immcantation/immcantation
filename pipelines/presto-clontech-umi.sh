#!/usr/bin/env bash
# BuildConsensus run parameters#!/usr/bin/env bash
# Script to run a pRESTO workflow on Takara SMARTer Human BCR IgG IgM H/K/L Profiling Kit
# https://www.takarabio.com/products/next-generation-sequencing/immune-profiling/human-repertoire/human-bcr-profiling-kit-for-illumina-sequencing
#
# Author:  Susanna Marquez, Jason Anthony Vander Heiden
# Date:    2021.09.11
#
# Arguments:
#   -1  Read 1 FASTQ sequence file (sequence beginning with the C-region or J-segment).
#   -2  Read 2 FASTQ sequence file (sequence beginning with the leader or V-segment).
#   -j  C-region reference sequences (reverse complemented).
#       Defaults to /usr/local/share/protocols/Universal/Human_IG_CRegion_RC.fasta
#   -r  V-segment reference file.
#       Defaults to /usr/local/share/igblast/fasta/imgt_human_ig_v.fasta
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
            "     Defaults to /usr/local/share/protocols/Universal/Human_IG_CRegion_RC.fasta."
    echo -e "  -r  V-segment reference file.\n" \
            "     Defaults to /usr/local/share/igblast/fasta/imgt_human_ig_v.fasta."
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
while getopts "1:2:j:r:y:n:o:x:p:h" OPT; do
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

# Check R1 reads
if [ -e ${R1_READS} ]; then
    R1_READS=$(realpath ${R1_READS})
else
    echo -e "File '${R1_READS}' not found." >&2
    exit 1
fi

# Check R2 reads
if [ -e ${R2_READS} ]; then
    R2_READS=$(realpath ${R2_READS})
else
    echo -e "File '${R2_READS}' not found." >&2
    exit 1
fi

# Check R1 primers
if ! ${C_PRIMERS_SET}; then
    C_PRIMERS="/usr/local/share/protocols/Universal/Human_IG_CRegion_RC.fasta"
elif [ -e ${C_PRIMERS} ]; then
    C_PRIMERS=$(realpath ${C_PRIMERS})
else
    echo -e "File '${C_PRIMERS}' not found." >&2
    exit 1
fi

# Check reference sequences
if ! ${VREF_SEQ_SET}; then
    VREF_SEQ="/usr/local/share/igblast/fasta/imgt_human_ig_v.fasta"
elif [ -e ${VREF_SEQ} ]; then
    VREF_SEQ=$(realpath ${VREF_SEQ})
else
    echo -e "File '${VREF_SEQ}' not found." >&2
    exit 1
fi

# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true
REPORT=true

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
MP_MAXERR=0.2
MP_MAXLEN=70
C_FIELD="C_CALL"

# AlignSets run parameters
MUSCLE_EXEC=muscle

# CollapseSeq run parameters
CS_KEEP=true
CS_MISS=0

# BuildConsensus run parameters
BC_QUAL=0
BC_MINCOUNT=1
BC_MAXERR=0.1
BC_PRCONS=0.6
BC_MAXGAP=0.5

# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="logs"
REPORTDIR="report"
PIPELINE_LOG="${LOGDIR}/pipeline-presto.log"
ERROR_LOG="${LOGDIR}/pipeline-presto.err"
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
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
FilterSeq.py quality -s $R1_READS -q $FS_QUAL --nproc $NPROC \
    --outname "${OUTNAME}-R1" --outdir . --log "${LOGDIR}/quality-R1.log" \
    >> $PIPELINE_LOG  2> $ERROR_LOG
FilterSeq.py quality -s $R2_READS -q $FS_QUAL --nproc $NPROC \
    --outname "${OUTNAME}-R2" --outdir . --log "${LOGDIR}/quality-R2.log" \
    >> $PIPELINE_LOG  2> $ERROR_LOG
MPR1_FILE="${OUTNAME}-R1_quality-pass.fastq"
MPR2_FILE="${OUTNAME}-R2_quality-pass.fastq"
check_error

# Identify primers and UMI in -2 reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers extract"
MaskPrimers.py extract -s ${MPR2_FILE} \
     --start 12 --len 7 --barcode --bf BARCODE --mode cut \
     --log "${LOGDIR}/primers-2.log" \
    --outname "${OUTNAME}-R2" --outdir .  >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Annotate -1 reads with internal C-region
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s ${MPR1_FILE} \
    -p ${C_PRIMERS} \
    --maxlen ${MP_MAXLEN} --maxerror ${MP_MAXERR} \
    --mode cut --skiprc --pf ${C_FIELD} \
    --log "${LOGDIR}/primers-1.log" \
    --outname "${OUTNAME}-R1" --nproc ${NPROC} \
    --outdir .  >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Transfer annotation
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 "${OUTNAME}-R2_primers-pass.fastq" \
    -2 "${OUTNAME}-R1_primers-pass.fastq" \
    --1f BARCODE --2f ${C_FIELD} --coord ${COORD}  >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Multiple align UID read groups
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AlignSets"
AlignSets.py muscle -s "${OUTNAME}-R1_primers-pass_pair-pass.fastq" --exec $MUSCLE_EXEC \
        --nproc ${NPROC} --log "${LOGDIR}/align-1.log" \
        --outname "${OUTNAME}-R1"  >> $PIPELINE_LOG 2> $ERROR_LOG
AlignSets.py muscle -s "${OUTNAME}-R2_primers-pass_pair-pass.fastq" --exec $MUSCLE_EXEC \
        --nproc ${NPROC} --log "${LOGDIR}/align-2.log" \
        --outname "${OUTNAME}-R2"  >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# UMI consensus
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "BuildConsensus"
BuildConsensus.py -s "${OUTNAME}-R2_align-pass.fastq" \
   --bf BARCODE --pf ${C_FIELD} --prcons ${BC_PRCONS} \
   -n ${BC_MINCOUNT} -q ${BC_QUAL} --maxerror ${BC_MAXERR} --maxgap ${BC_MAXGAP}  \
   --nproc ${NPROC} --log "${LOGDIR}/consensus-2.log" \
   --outdir . --outname "${OUTNAME}-R2" >> $PIPELINE_LOG 2> $ERROR_LOG
BuildConsensus.py -s "${OUTNAME}-R1_align-pass.fastq" \
   --bf BARCODE --pf ${C_FIELD} --prcons ${BC_PRCONS} \
   -n ${BC_MINCOUNT} -q ${BC_QUAL} --maxerror ${BC_MAXERR} --maxgap ${BC_MAXGAP}  \
   --nproc ${NPROC} --log "${LOGDIR}/consensus-1.log"  \
   --outdir . --outname "${OUTNAME}-R1"  >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Syncronize reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 "${OUTNAME}-R1_consensus-pass.fastq" \
    -2 "${OUTNAME}-R2_consensus-pass.fastq" \
    --coord presto >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Assemble pairs
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs sequential"
AssemblePairs.py sequential -1 "${OUTNAME}-R1_consensus-pass_pair-pass.fastq" \
    -2 "${OUTNAME}-R2_consensus-pass_pair-pass.fastq" \
    -r ${VREF_SEQ} \
    --coord presto --rc tail --1f CONSCOUNT --2f PRCONS CONSCOUNT \
    --minlen ${AP_MINLEN} --maxerror ${AP_MAXERR} --alpha $AP_ALPHA --scanrev \
    --minident ${AP_MINIDENT} --evalue ${AP_EVALUE} --maxhits ${AP_MAXHITS} --aligner blastn \
    --nproc  $NPROC --log "${LOGDIR}/assemble.log" \
    --outname "${OUTNAME}" >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Mask low quality positions
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq maskqual"
FilterSeq.py maskqual -s "${OUTNAME}_assemble-pass.fastq" -q ${FS_MASK} --nproc ${NPROC} \
        --outname "${OUTNAME}" --log "${LOGDIR}/maskqual.log" >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Rewrite header with minimum of CONSCOUNT
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders collapse"
ParseHeaders.py collapse -s "${OUTNAME}_maskqual-pass.fastq" -f CONSCOUNT --act min \
    --outname "${OUTNAME}-final-prcons" >> $PIPELINE_LOG 2> $ERROR_LOG
mv "${OUTNAME}-final-prcons_reheader.fastq" "${OUTNAME}-final-prcons_total.fastq"
check_error

# Rename PRCONS to C_CALL
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders rename"
ParseHeaders.py rename -s "${OUTNAME}-final-prcons_total.fastq" \
    -f PRCONS -k C_CALL \
    --outname "${OUTNAME}-final_total" >> $PIPELINE_LOG 2> $ERROR_LOG
mv "${OUTNAME}-final_total_reheader.fastq" "${OUTNAME}-final_total.fastq"
check error

# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
CollapseSeq.py -s "${OUTNAME}-final_total.fastq" -n ${CS_MISS} \
    --uf C_CALL --cf CONSCOUNT --act sum --inner \
    --keepmiss --outname "${OUTNAME}-final" >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Subset to sequences seen at least twice
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq"
SplitSeq.py group -s "${OUTNAME}-final_collapse-unique.fastq" \
    -f CONSCOUNT --num 2 >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Create table of final repertoire
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
ParseHeaders.py table -s "${OUTNAME}-final_total.fastq" -f ID ${C_FIELD} CONSCOUNT \
    --outname "final-total" --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
ParseHeaders.py table -s "${OUTNAME}-final_collapse-unique.fastq" -f ID ${C_FIELD} CONSCOUNT DUPCOUNT \
    --outname "final-unique" --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
ParseHeaders.py table -s "${OUTNAME}-final_collapse-unique_atleast-2.fastq" -f ID ${C_FIELD} CONSCOUNT DUPCOUNT \
    --outname "final-unique-atleast2" --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
check_error


# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l "logs/quality-R1.log" "logs/quality-R2.log" -f ID QUALITY \
        --outdir ${LOGDIR} > /dev/null     
ParseLog.py -l "${LOGDIR}/primers-1.log" -f ID PRIMER ERROR PRSTART \
        --outdir "${LOGDIR}" > /dev/null 
ParseLog.py -l "${LOGDIR}/align-1.log" "${LOGDIR}/align-2.log" -f BARCODE SEQCOUNT \
        --outdir "${LOGDIR}" > /dev/null 
ParseLog.py -l "${LOGDIR}/consensus-1.log" "${LOGDIR}/consensus-2.log" \
    -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ ERROR \
        --outdir "${LOGDIR}" > /dev/null 
ParseLog.py -l "${LOGDIR}/assemble.log" \
    -f ID REFID LENGTH OVERLAP GAP ERROR PVALUE EVALUE1 EVALUE2 IDENTITY FIELDS1 FIELDS2 \
        --outdir "${LOGDIR}" > /dev/null 
ParseLog.py -l "${LOGDIR}/maskqual.log" -f ID MASKED \
        --outdir "${LOGDIR}" > /dev/null 

wait
check_error

# Generate pRESTO report
if $REPORT; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Generating report"
    REPORT_SCRIPT="buildReport(\"${LOGDIR}\", sample=\"${OUTNAME}\", output_dir=\"${REPORTDIR}\", template=\"Clontech-UMI\", config=NULL, quiet=FALSE)"
    Rscript -e "library(prestor); ${REPORT_SCRIPT}" > ${REPORTDIR}/report.out 2> ${REPORTDIR}/report.err
fi

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
cd ..
