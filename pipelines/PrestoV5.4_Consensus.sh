#!/usr/bin/env bash
# Super script to run consensus building steps of the pRESTO 0.5.4 pipeline on AbVitro AbSeq data
# 
# Author:  Jason Anthony Vander Heiden, Gur Yaari, Namita Gupta
# Date:    2017.10.24
# 
# Arguments:
#   -s  FASTQ sequence file.
#   -c  C-region FASTA sequences for the C-region internal to the primer.
#       If unspecified internal C-region alignment is not performed.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the read 1 filename.
#   -o  Output directory.
#       Defaults to the sample name.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -s  FASTQ sequence file."
    echo -e "  -c  C-region FASTA sequences for the C-region internal to the primer.\n" \
            "     If unspecified internal C-region alignment is not performed."
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the read 1 filename."
    echo -e "  -o  Output directory.\n" \
            "     Defaults to the sample name."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -h  This message."
}

# Argument validation variables
READS_SET=false
CREGION_SEQ_SET=false
OUTNAME_SET=false
OUTDIR_SET=false
NPROC_SET=false

# Get commandline arguments
while getopts "s:c::n:o:p:h" OPT; do
    case "$OPT" in
    s)  READS=${OPTARG}
        READS_SET=true
        ;;
    c)  CREGION_SEQ=${OPTARG}
        CREGION_SEQ_SET=true
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
if ! ${READS_SET}; then
    echo -e "You must specify a read file using the -s option." >&2
    exit 1
elif [ -e ${READS} ]; then
    READS=$(readlink -f ${READS})
else
    echo -e "File ${READS} not found." >&2
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

# Check for C-region file
if ! ${CREGION_SEQ_SET}; then
    ALIGN_CREGION=false
elif [ -e ${CREGION_SEQ} ]; then
    ALIGN_CREGION=true
    CREGION_SEQ=$(readlink -f ${CREGION_SEQ})
else
    echo -e "File ${CREGION_SEQ} not found." >&2
    exit 1
fi

# Define pipeline steps
ZIP_FILES=true
DELETE_FILES=true
ALIGN_SETS=false
MASK_LOWQUAL=false

# FilterSeq run parameters
FS_MASK=30

# MaskPrimers run parameters
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
if $BC_ERR_FLAG; then
    BC_ERROR="--maxerror ${BC_MAXERR}"
else
    BC_ERROR=""
fi

# CollapseSeq run parameters
CS_KEEP=false
CS_MISS=20
if $CS_KEEP; then
    CS_KEEPMISS="--keepmiss"
else
    CS_KEEPMISS=""
fi

# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="logs"
PIPELINE_LOG="${LOGDIR}/pipeline-consensus.log"
ERROR_LOG="${LOGDIR}/pipeline-consensus.err"
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

# Multiple align UID read groups
if $ALIGN_SETS; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AlignSets muscle"
	AlignSets.py muscle -s $READS --exec $MUSCLE_EXEC \
	    --nproc $NPROC --log "${LOGDIR}/align.log" --outname "${OUTNAME}" --outdir . \
	    >> $PIPELINE_LOG 2> $ERROR_LOG
	    
	BC_FILE="${OUTNAME}_align-pass.fastq"
	check_error
else
	BC_FILE=$READS
fi

if $BC_PRCONS_FLAG; then
	printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders expand"
    ParseHeaders.py expand -s $BC_FILE -f PRIMER \
        --outname "${OUTNAME}" --outdir . >> $PIPELINE_LOG 2> $ERROR_LOG

	printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders rename"
    ParseHeaders.py rename -s "${OUTNAME}_reheader.fastq" -f PRIMER1 PRIMER2 \
        -k FPRIMER RPRIMER >> $PIPELINE_LOG 2> $ERROR_LOG

    BC_FILE="${OUTNAME}_reheader_reheader.fastq"
    check_error
fi

# Build UID consensus sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "BuildConsensus"
if $BC_PRCONS_FLAG; then
    BuildConsensus.py -s $BC_FILE --bf BARCODE --pf RPRIMER --prcons $BC_PRCONS \
        -n $BC_MINCOUNT -q $BC_QUAL --maxgap $BC_MAXGAP $BC_ERROR \
        --nproc $NPROC --log "${LOGDIR}/consensus.log" \
        --outname "${OUTNAME}" >> $PIPELINE_LOG 2> $ERROR_LOG
else
    BuildConsensus.py -s $BC_FILE --bf BARCODE \
        -n $BC_MINCOUNT -q $BC_QUAL --maxgap $BC_MAXGAP $BC_ERROR \
        --nproc $NPROC --log "${LOGDIR}/consensus.log" \
        --outname "${OUTNAME}" --outdir . >> $PIPELINE_LOG 2> $ERROR_LOG
fi
check_error

# Mask low quality positions
if $MASK_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq maskqual"
    FilterSeq.py maskqual -s "${OUTNAME}_consensus-pass.fastq" -q $FS_MASK --nproc $NPROC \
        --outname "${OUTNAME}" --log "${LOGDIR}/maskqual.log" \
        >> $PIPELINE_LOG 2> $ERROR_LOG

    PH_FILE="${OUTNAME}_maskqual-pass.fastq"
    check_error
else
    PH_FILE="${OUTNAME}_consensus-pass.fastq"
fi

if $ALIGN_CREGION; then
    # Annotate with internal C-region
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
    MaskPrimers.py align -s $PH_FILE -p $CREGION_SEQ \
        --maxlen $CREGION_MAXLEN --maxerror $CREGION_MAXERR \
        --mode tag --revpr --skiprc \
        --log "${LOGDIR}/cregion.log" --outname "${OUTNAME}" --nproc $NPROC \
        >> $PIPELINE_LOG 2> $ERROR_LOG

    # Renamer primer field
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders rename"
    ParseHeaders.py rename -s "${OUTNAME}_primers-pass.fastq" -f PRIMER -k CREGION \
        --outname "${OUTNAME}" > /dev/null 2> $ERROR_LOG
    PH_FILE="${OUTNAME}_reheader.fastq"

    CREGION_FIELD="CREGION"
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
CollapseSeq.py -s "${OUTNAME}-final_total.fastq" -n $CS_MISS \
    --uf PRCONS $CREGION_FIELD --cf CONSCOUNT --act sum --inner \
    ${CS_KEEPMISS} --outname "${OUTNAME}-final" >> $PIPELINE_LOG 2> $ERROR_LOG
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
if $ALIGN_SETS; then
    ParseLog.py -l "${LOGDIR}/align.log" -f BARCODE SEQCOUNT \
        --outdir ${LOGDIR} > /dev/null  2> $ERROR_LOG &
fi
ParseLog.py -l "${LOGDIR}/consensus.log" \
    -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ ERROR \
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

# Zip or delete intermediate and log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Compressing files"
LOG_FILES=$(ls ${LOGDIR}/*.log | grep -v "pipeline")
FILTER_FILES="$(basename ${READS})"
FILTER_FILES+="\|final_total.fastq\|final_collapse-unique.fastq\|final_collapse-unique_atleast-2.fastq"
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
