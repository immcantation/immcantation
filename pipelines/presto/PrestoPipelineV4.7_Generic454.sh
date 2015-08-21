#!/bin/bash
# Super script to run the pRESTO 0.4.7 pipeline on 454 data
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2015.08.20
# 
# Required Arguments:
#   $1 = read file
#   $2 = MID file
#   $3 = forward (V-region) primer file
#   $4 = reverse (C-region) primer file
#   $5 = output directory
#   $6 = output file prefix
#   $7 = number of subprocesses for multiprocessing tools

# Capture command line parameters
READ_FILE=$(readlink -f $1)
MID_PRIMERS=$(readlink -f $2)
FWD_PRIMERS=$(readlink -f $3)
REV_PRIMERS=$(readlink -f $4)
OUTDIR=$5
OUTNAME=$6
NPROC=$7

# Define pipeline steps
ZIP_FILES=true
FILTER_LENGTH=true
FILTER_LOWQUAL=true
MASK_LOWQUAL=false

# FilterSeq run parameters
FS_LENGTH=300
FS_QUAL=20
FS_MASK=30

# MaskPrimers run parameters
MP_MID_START=0
MP_MID_MODE="cut"
MP_FWD_MODE="mask"
MP_REV_MODE="cut"
MP_MID_MAXERR=0.2
MP_FWD_MAXERR=0.3
MP_REV_MAXERR=0.3
MP_FWD_MAXLEN=50
MP_REV_MAXLEN=50

# CollapseSeq run parameters
CS_KEEP=false
CS_MISS=20

# Define log files
PIPELINE_LOG="Pipeline.log"
ERROR_LOG="Pipeline.err"

# Make output directory and empty log files
mkdir -p $OUTDIR; cd $OUTDIR
echo '' > $PIPELINE_LOG
echo '' > $ERROR_LOG

# Start
echo "DIRECTORY: ${OUTDIR}"
echo "VERSIONS:"
echo "  $(AlignSets.py --version 2>&1)"
echo "  $(AssemblePairs.py --version 2>&1)"
echo "  $(BuildConsensus.py --version 2>&1)"
echo "  $(ClusterSets.py --version 2>&1)"
echo "  $(CollapseSeq.py --version 2>&1)"
echo "  $(ConvertHeaders.py --version 2>&1)"
echo "  $(FilterSeq.py --version 2>&1)"
echo "  $(MaskPrimers.py --version 2>&1)"
echo "  $(PairSeq.py --version 2>&1)"
echo "  $(ParseHeaders.py --version 2>&1)"
echo "  $(ParseLog.py --version 2>&1)"
echo "  $(SplitSeq.py --version 2>&1)"
echo -e "\nSTART"
STEP=0

# Remove short reads
if $FILTER_LENGTH; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq length"
    FilterSeq.py length -s $READ_FILE -n $FS_LENGTH --nproc $NPROC \
        --outname "${OUTNAME}" --outdir . --log LengthLog.log \
        >> $PIPELINE_LOG  2> $ERROR_LOG
    MP_FILE="${OUTNAME}_length-pass.fastq"
else
    MP_FILE=$READ_FILE
fi

# Remove low quality reads
if $FILTER_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
    FilterSeq.py quality -s $MP_FILE -q $FS_QUAL --nproc $NPROC \
        --outname "${OUTNAME}" --outdir . --log QualityLog.log \
        >> $PIPELINE_LOG  2> $ERROR_LOG
    MP_FILE="${OUTNAME}_quality-pass.fastq"
fi

# Identify and remove MIDs
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
MaskPrimers.py score -s $MP_FILE -p $MID_PRIMERS --mode $MP_MID_MODE \
    --start $MP_MID_START --maxerror $MP_MID_MAXERR --nproc $NPROC \
    --outname "${OUTNAME}" --outdir . --log PrimerMIDLog.log \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Remove forward (V-region) primers
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s "${OUTNAME}_primers-pass.fastq" -p $FWD_PRIMERS \
    --mode $MP_FWD_MODE --maxlen $MP_FWD_MAXLEN --maxerror $MP_FWD_MAXERR \
    --nproc $NPROC --log PrimerForwardLog.log \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Remove reverse (C-region) primers
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s "${OUTNAME}_primers-pass_primers-pass.fastq" -p $REV_PRIMERS \
    --mode $MP_REV_MODE --maxlen $MP_REV_MAXLEN --maxerror $MP_REV_MAXERR \
    --revpr --skiprc --nproc $NPROC --log PrimerReverseLog.log \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Expand primer field
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders expand"
ParseHeaders.py expand -s "${OUTNAME}_primers-pass_primers-pass_primers-pass.fastq" \
    -f PRIMER --outname "${OUTNAME}-FIN" \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Rename primer fields
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders expand"
ParseHeaders.py rename -s "${OUTNAME}-FIN_reheader.fastq" \
    -f PRIMER1 PRIMER2 PRIMER3 -k MID VPRIMER CPRIMER \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Mask low quality positions
if $MASK_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq maskqual"
    FilterSeq.py maskqual -s "${OUTNAME}-FIN_reheader_reheader.fastq" -q $FS_MASK \
        --nproc $NPROC --log MaskQualLog.log \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    SS_FILE="${OUTNAME}-FIN_reheader_reheader_maskqual-pass.fastq"
else
    SS_FILE="${OUTNAME}-FIN_reheader_reheader.fastq"
fi

# Split file by MID
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s $SS_FILE -f MID --outname "${OUTNAME}-FIN" \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
if $CS_KEEP; then
    CollapseSeq.py -s "${OUTNAME}-FIN_MID-"* -n $CS_MISS --uf CPRIMER --inner --keepmiss \
    >> $PIPELINE_LOG 2> $ERROR_LOG
else
    CollapseSeq.py -s "${OUTNAME}-FIN_MID-"* -n $CS_MISS --uf CPRIMER --inner \
    >> $PIPELINE_LOG 2> $ERROR_LOG
fi

# Filter to atleast 2 reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s *_collapse-unique.fastq -f DUPCOUNT --num 2 \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
if $FILTER_LENGTH; then
    ParseLog.py -l LengthLog.log -f ID LENGTH \
    > /dev/null 2> $ERROR_LOG &
fi
if $FILTER_LOWQUAL; then
    ParseLog.py -l QualityLog.log -f ID QUALITY \
    > /dev/null 2> $ERROR_LOG &
fi
if $MASK_LOWQUAL; then
    ParseLog.py -l MaskQualLog.log -f ID MASKED \
    > /dev/null  2> $ERROR_LOG &
fi
ParseLog.py -l PrimerMIDLog.log PrimerForwardLog.log PrimerReverseLog.log \
    -f ID PRIMER ERROR \
    > /dev/null  2> $ERROR_LOG &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    LOG_FILES_ZIP=$(ls *Log.log)
    tar -zcf LogFiles.tar $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP

    TEMP_FILES_ZIP=$(ls *.fastq | grep -vP "collapse-unique.fastq|collapse-unique_atleast-2.fastq")
    tar -zcf TempFiles.tar $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
fi

# End
printf "DONE\n\n"
cd ../

