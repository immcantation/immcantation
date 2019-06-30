#!/usr/bin/env bash
# Super script to run the pRESTO pipeline on AbVitro/NEB AbSeq data
# 
# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2018.09.15
# 
# Arguments:
#   -1  Read 1 FASTQ sequence file (sequence beginning with the C-region or J-segment).
#   -2  Read 2 FASTQ sequence file (sequence beginning with the leader or V-segment).
#   -j  Read 1 FASTA primer sequences (C-region or J-segment).
#       Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R1_Human_IG_Primers.fasta
#   -v  Read 2 FASTA primer sequences (template switch or V-segment).
#       Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R2_TS.fasta.
#   -c  C-region FASTA sequences for the C-region internal to the primer.
#       If unspecified internal C-region alignment is not performed.
#   -r  V-segment reference file.
#       Defaults to /usr/local/share/germlines/igblast/fasta/imgt_human_ig_v.fasta
#   -y  YAML file providing description fields for report generation.
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
            "     Sequence beginning with the C-region or J-segment)."
    echo -e "  -2  Read 2 FASTQ sequence file.\n" \
            "     Sequence beginning with the leader or V-segment)."
    echo -e "  -j  Read 1 FASTA primer sequences.\n" \
            "     Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R1_Human_IG_Primers.fasta."
    echo -e "  -v  Read 2 FASTA primer or template switch sequences.\n" \
            "     Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R2_TS.fasta."
    echo -e "  -c  C-region FASTA sequences for the C-region internal to the primer.\n" \
            "     If unspecified internal C-region alignment is not performed."
    echo -e "  -r  V-segment reference file.\n" \
            "     Defaults to /usr/local/share/igblast/fasta/imgt_human_ig_v.fasta."
    echo -e "  -y  YAML file providing description fields for report generation."
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
R1_PRIMERS_SET=false
R2_PRIMERS_SET=false
CREGION_SEQ_SET=false
VREF_SEQ_SET=false
YAML_SET=FALSE
OUTNAME_SET=false
OUTDIR_SET=false
NPROC_SET=false
COORD_SET=false

# Get commandline arguments
while getopts "1:2:j:v:c:r:y:n:o:x:p:h" OPT; do
    case "$OPT" in
    1)  R1_READS=$OPTARG
        R1_READS_SET=true
        ;;
    2)  R2_READS=$OPTARG
        R2_READS_SET=true
        ;;
    j)  R1_PRIMERS=$OPTARG
        R1_PRIMERS_SET=true
        ;;
    v)  R2_PRIMERS=$OPTARG
        R2_PRIMERS_SET=true
        ;;
    c)  CREGION_SEQ=$OPTARG
        CREGION_SEQ_SET=true
        ;;
    r)  VREF_SEQ=$OPTARG
        VREF_SEQ_SET=true
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
if ! ${R1_READS_SET} || ! ${R2_READS_SET}; then
    echo -e "You must specify both read files using the -1 and -2 options." >&2
    exit 1
fi

if ! ${YAML_SET}; then
    echo -e "You must specify the description file in YAML format using the -y option." >&2
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
    PARENTDIR=$(dirname $(readlink -f ${OUTDIR}))
    if ! [ -w ${PARENTDIR} ]; then
        echo -e "Parent directory '${PARENTDIR}' of new output directory '${OUTDIR}' is not writable." >&2
        exit 1
    fi
fi

# Check R1 reads
if [ -e ${R1_READS} ]; then
    R1_READS=$(readlink -f ${R1_READS})
else
    echo -e "File '${R1_READS}' not found." >&2
    exit 1
fi

# Check R2 reads
if [ -e ${R2_READS} ]; then
    R2_READS=$(readlink -f ${R2_READS})
else
    echo -e "File '${R2_READS}' not found." >&2
    exit 1
fi

# Check R1 primers
if ! ${R1_PRIMERS_SET}; then
    R1_PRIMERS="/usr/local/share/protocols/AbSeq/AbSeq_R1_Human_IG_Primers.fasta"
elif [ -e ${R1_PRIMERS} ]; then
    R1_PRIMERS=$(readlink -f ${R1_PRIMERS})
else
    echo -e "File '${R1_PRIMERS}' not found." >&2
    exit 1
fi

# Check R2 primers
if ! ${R2_PRIMERS_SET}; then
    R2_PRIMERS="/usr/local/share/protocols/AbSeq/AbSeq_R2_TS.fasta"
elif [ -e ${R2_PRIMERS} ]; then
    R2_PRIMERS=$(readlink -f ${R2_PRIMERS})
else
    echo -e "File '${R2_PRIMERS}' not found." >&2
    exit 1
fi

# Check reference sequences
if ! ${VREF_SEQ_SET}; then
    VREF_SEQ="/usr/local/share/igblast/fasta/imgt_human_ig_v.fasta"
elif [ -e ${VREF_SEQ} ]; then
    VREF_SEQ=$(readlink -f ${VREF_SEQ})
else
    echo -e "File '${VREF_SEQ}' not found." >&2
    exit 1
fi

# Check for C-region file
if ! ${CREGION_SEQ_SET}; then
    ALIGN_CREGION=false
elif [ -e ${CREGION_SEQ} ]; then
    ALIGN_CREGION=true
    CREGION_SEQ=$(readlink -f ${CREGION_SEQ})
else
    echo -e "File '${CREGION_SEQ}' not found." >&2
    exit 1
fi

# Check report yaml file
if [ -e ${YAML} ]; then
    YAML=$(readlink -f ${YAML})
else
    echo -e "File '${YAML}' not found." >&2
    exit 1
fi

# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true
FILTER_LOWQUAL=true
ALIGN_SETS=false
MASK_LOWQUAL=false
REPORT=true

# FilterSeq run parameters
FS_QUAL=20
FS_MASK=30

# MaskPrimers run parameters
MP_UIDLEN=17
MP_R1_MAXERR=0.2
MP_R2_MAXERR=0.5
CREGION_MAXLEN=100
CREGION_MAXERR=0.3

# AlignSets run parameters
MUSCLE_EXEC=muscle

# BuildConsensus run parameters
BC_PRCONS_FLAG=true
BC_ERR_FLAG=true
BC_QUAL=0
BC_MINCOUNT=1
BC_MAXERR=0.1
BC_PRCONS=0.6
BC_MAXGAP=0.5

# AssemblePairs-sequential run parameters
AP_MAXERR=0.3
AP_MINLEN=8
AP_ALPHA=1e-5
AP_MINIDENT=0.5
AP_EVALUE=1e-5
AP_MAXHITS=100

# CollapseSeq run parameters
CS_KEEP=true
CS_MISS=0

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
if $FILTER_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
    FilterSeq.py quality -s $R1_READS -q $FS_QUAL --nproc $NPROC \
        --outname "${OUTNAME}-R1" --outdir . --log "${LOGDIR}/quality-1.log" \
        >> $PIPELINE_LOG  2> $ERROR_LOG
    FilterSeq.py quality -s $R2_READS -q $FS_QUAL --nproc $NPROC \
        --outname "${OUTNAME}-R2" --outdir . --log "${LOGDIR}/quality-2.log"  \
        >> $PIPELINE_LOG  2> $ERROR_LOG
    MPR1_FILE="${OUTNAME}-R1_quality-pass.fastq"
    MPR2_FILE="${OUTNAME}-R2_quality-pass.fastq"
    check_error
else
    MPR1_FILE=$R1_READS
    MPR2_FILE=$R2_READS
fi


# Identify primers and UID 
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
MaskPrimers.py score -s $MPR1_FILE -p $R1_PRIMERS --mode cut \
    --start 0 --maxerror $MP_R1_MAXERR --nproc $NPROC \
    --log "${LOGDIR}/primers-1.log" --outname "${OUTNAME}-R1" --outdir . \
    >> $PIPELINE_LOG 2> $ERROR_LOG
MaskPrimers.py score -s $MPR2_FILE -p $R2_PRIMERS --mode cut \
    --start $MP_UIDLEN --barcode --maxerror $MP_R2_MAXERR --nproc $NPROC \
    --log "${LOGDIR}/primers-2.log" --outname "${OUTNAME}-R2" --outdir . \
    >> $PIPELINE_LOG 2> $ERROR_LOG
check_error


# Assign UIDs to read 1 sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 "${OUTNAME}-R1_primers-pass.fastq" -2 "${OUTNAME}-R2_primers-pass.fastq" \
    --2f BARCODE --coord $COORD >> $PIPELINE_LOG 2> $ERROR_LOG
check_error


# Multiple align UID read groups
if $ALIGN_SETS; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AlignSets muscle"
	AlignSets.py muscle -s "${OUTNAME}-R1_primers-pass_pair-pass.fastq" --exec $MUSCLE_EXEC \
	    --nproc $NPROC --log "${LOGDIR}/align-1.log" --outname "${OUTNAME}-R1" \
	    >> $PIPELINE_LOG 2> $ERROR_LOG
	AlignSets.py muscle -s "${OUTNAME}-R2_primers-pass_pair-pass.fastq" --exec $MUSCLE_EXEC \
	    --nproc $NPROC --log "${LOGDIR}/align-2.log" --outname "${OUTNAME}-R2" \
	    >> $PIPELINE_LOG 2> $ERROR_LOG
	BCR1_FILE="${OUTNAME}-R1_align-pass.fastq"
	BCR2_FILE="${OUTNAME}-R2_align-pass.fastq"
	check_error
else
	BCR1_FILE="${OUTNAME}-R1_primers-pass_pair-pass.fastq"
	BCR2_FILE="${OUTNAME}-R2_primers-pass_pair-pass.fastq"
fi


# Build UID consensus sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "BuildConsensus"
if $BC_ERR_FLAG; then
    if $BC_PRCONS_FLAG; then
        BuildConsensus.py -s $BCR1_FILE --bf BARCODE --pf PRIMER --prcons $BC_PRCONS \
            -n $BC_MINCOUNT -q $BC_QUAL --maxerror $BC_MAXERR --maxgap $BC_MAXGAP \
            --nproc $NPROC --log "${LOGDIR}/consensus-1.log" \
            --outname "${OUTNAME}-R1" >> $PIPELINE_LOG 2> $ERROR_LOG
    else
        BuildConsensus.py -s $BCR1_FILE --bf BARCODE --pf PRIMER \
            -n $BC_MINCOUNT -q $BC_QUAL --maxerror $BC_MAXERR --maxgap $BC_MAXGAP \
            --nproc $NPROC --log "${LOGDIR}/consensus-1.log" \
            --outname "${OUTNAME}-R1" >> $PIPELINE_LOG 2> $ERROR_LOG
    fi

	BuildConsensus.py -s $BCR2_FILE --bf BARCODE --pf PRIMER \
	    -n $BC_MINCOUNT -q $BC_QUAL --maxerror $BC_MAXERR --maxgap $BC_MAXGAP \
	    --nproc $NPROC --log "${LOGDIR}/consensus-2.log" \
	    --outname "${OUTNAME}-R2" >> $PIPELINE_LOG 2> $ERROR_LOG
else
    if $BC_PRCONS_FLAG; then
        BuildConsensus.py -s $BCR1_FILE --bf BARCODE --pf PRIMER --prcons $BC_PRCONS \
            -n $BC_MINCOUNT -q $BC_QUAL --maxgap $BC_MAXGAP \
            --nproc $NPROC --log "${LOGDIR}/consensus-1.log" \
            --outname "${OUTNAME}-R1" >> $PIPELINE_LOG 2> $ERROR_LOG
    else
        BuildConsensus.py -s $BCR1_FILE --bf BARCODE --pf PRIMER \
            -n $BC_MINCOUNT -q $BC_QUAL --maxgap $BC_MAXGAP \
            --nproc $NPROC --log "${LOGDIR}/consensus-1.log" \
            --outname "${OUTNAME}-R1" >> $PIPELINE_LOG 2> $ERROR_LOG
    fi

	BuildConsensus.py -s $BCR2_FILE --bf BARCODE --pf PRIMER \
    	-n $BC_MINCOUNT -q $BC_QUAL --maxgap $BC_MAXGAP \
    	--nproc $NPROC --log "${LOGDIR}/consensus-2.log" \
    	--outname "${OUTNAME}-R2" >> $PIPELINE_LOG 2> $ERROR_LOG
fi
check_error


# Syncronize read files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 "${OUTNAME}-R1_consensus-pass.fastq" -2 "${OUTNAME}-R2_consensus-pass.fastq" \
    --coord presto >> $PIPELINE_LOG 2> $ERROR_LOG
check_error


# Assemble paired ends via mate-pair alignment
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs sequential"
if $BC_PRCONS_FLAG; then
    PRFIELD="PRCONS"
else
    PRFIELD="PRIMER"
fi

AssemblePairs.py sequential -1 "${OUTNAME}-R2_consensus-pass_pair-pass.fastq" \
    -2 "${OUTNAME}-R1_consensus-pass_pair-pass.fastq" -r $VREF_SEQ \
    --coord presto --rc tail --1f CONSCOUNT --2f $PRFIELD CONSCOUNT \
    --minlen $AP_MINLEN --maxerror $AP_MAXERR --alpha $AP_ALPHA --scanrev \
    --minident $AP_MINIDENT --evalue $AP_EVALUE --maxhits $AP_MAXHITS --aligner blastn \
    --nproc $NPROC --log "${LOGDIR}/assemble.log" \
    --outname "${OUTNAME}" >> $PIPELINE_LOG 2> $ERROR_LOG
PH_FILE="${OUTNAME}_assemble-pass.fastq"
check_error


# Mask low quality positions
if $MASK_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq maskqual"
    FilterSeq.py maskqual -s $PH_FILE -q $FS_MASK --nproc $NPROC \
        --outname "${OUTNAME}-MQ" --log "${LOGDIR}/maskqual.log" \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    PH_FILE="${OUTNAME}-MQ_maskqual-pass.fastq"
    check_error
fi


if $ALIGN_CREGION; then
    # Annotate with internal C-region
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
    CREGION_FIELD="CREGION"
    MaskPrimers.py align -s $PH_FILE -p $CREGION_SEQ \
        --maxlen $CREGION_MAXLEN --maxerror $CREGION_MAXERR \
        --mode tag --revpr --skiprc --pf $CREGION_FIELD \
        --log "${LOGDIR}/cregion.log" --outname "${OUTNAME}-CR" --nproc $NPROC \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    PH_FILE="${OUTNAME}-CR_primers-pass.fastq"
    check_error
else
    CREGION_FIELD=""
fi


# Rewrite header with minimum of CONSCOUNT
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders collapse"
ParseHeaders.py collapse -s $PH_FILE -f CONSCOUNT --act min \
    --outname "${OUTNAME}-final" > /dev/null 2> $ERROR_LOG
mv "${OUTNAME}-final_reheader.fastq" "${OUTNAME}-final_total.fastq"
check_error


# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
if $CS_KEEP; then
    CollapseSeq.py -s "${OUTNAME}-final_total.fastq" -n $CS_MISS \
    --uf PRCONS $CREGION_FIELD --cf CONSCOUNT --act sum --inner \
    --keepmiss --outname "${OUTNAME}-final" >> $PIPELINE_LOG 2> $ERROR_LOG
else
    CollapseSeq.py -s "${OUTNAME}-final_total.fastq" -n $CS_MISS \
    --uf PRCONS $CREGION_FIELD --cf CONSCOUNT --act sum --inner \
    --outname "${OUTNAME}-final" >> $PIPELINE_LOG 2> $ERROR_LOG
fi
check_error


# Filter to sequences with at least 2 supporting sources
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s "${OUTNAME}-final_collapse-unique.fastq" -f CONSCOUNT --num 2 \
    >> $PIPELINE_LOG 2> $ERROR_LOG
check_error


# Create table of final repertoire
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
ParseHeaders.py table -s "${OUTNAME}-final_total.fastq" \
    -f ID PRCONS $CREGION_FIELD CONSCOUNT --outname "final-total" \
    --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
ParseHeaders.py table -s "${OUTNAME}-final_collapse-unique.fastq" \
    -f ID PRCONS $CREGION_FIELD CONSCOUNT DUPCOUNT --outname "final-unique" \
    --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
ParseHeaders.py table -s "${OUTNAME}-final_collapse-unique_atleast-2.fastq" \
    -f ID PRCONS $CREGION_FIELD CONSCOUNT DUPCOUNT --outname "final-unique-atleast2" \
    --outdir ${LOGDIR} >> $PIPELINE_LOG 2> $ERROR_LOG
check_error


# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
if $FILTER_LOWQUAL; then
    ParseLog.py -l "${LOGDIR}/quality-1.log" "${LOGDIR}/quality-2.log" -f ID QUALITY \
        --outdir ${LOGDIR} > /dev/null &
fi
ParseLog.py -l "${LOGDIR}/primers-1.log" "${LOGDIR}/primers-2.log" -f ID BARCODE PRIMER ERROR \
    --outdir ${LOGDIR} > /dev/null  2> $ERROR_LOG &
ParseLog.py -l "${LOGDIR}/consensus-1.log" "${LOGDIR}/consensus-2.log" \
    -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ ERROR \
    --outdir ${LOGDIR} > /dev/null  2> $ERROR_LOG &
ParseLog.py -l "${LOGDIR}/assemble.log" \
    -f ID REFID LENGTH OVERLAP GAP ERROR PVALUE EVALUE1 EVALUE2 IDENTITY FIELDS1 FIELDS2 \
    --outdir ${LOGDIR} > /dev/null  2> $ERROR_LOG &
if $MASK_LOWQUAL; then
    ParseLog.py -l "${LOGDIR}/maskqual.log" -f ID MASKED \
        --outdir ${LOGDIR} > /dev/null  2> $ERROR_LOG &
fi
if $ALIGN_CREGION; then
    ParseLog.py -l "${LOGDIR}/cregion.log" -f ID PRIMER ERROR \
        --outdir ${LOGDIR} > /dev/null  2> $ERROR_LOG &
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
FILTER_FILES="$(basename ${R1_READS})\|$(basename ${R2_READS})\|$(basename ${R1_PRIMERS})\|$(basename ${R2_PRIMERS})"
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

