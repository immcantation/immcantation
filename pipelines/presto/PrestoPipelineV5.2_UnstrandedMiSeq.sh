#!/bin/bash
# Super script to run the pRESTO 0.5.2 pipeline on unstranded MiSeq data without UIDs
#
# Author:  Jason Anthony Vander Heiden
# Date:    2015.03.10
#
# Required Arguments:
#   $1 = read 1 file (V-region start sequence)
#   $2 = read 2 file (C-region start sequence)
#   $3 = read 1 primer file (V-region)
#   $4 = read 2 primer file (C-region)
#   $5 = output directory
#   $6 = output file prefix
#   $7 = number of subprocesses for multiprocessing tools

# Capture command line parameters
R1_FILE=$(readlink -f $1)
R2_FILE=$(readlink -f $2)
R1_PRIMERS=$(readlink -f $3)
R2_PRIMERS=$(readlink -f $4)
OUTDIR=$5
OUTNAME=$6
NPROC=$7

# Define pipeline steps
ZIP_FILES=true
REFERENCE_ASSEMBLY=false
FILTER_LOWQUAL=true
MASK_LOWQUAL=false

# FilterSeq run parameters
FS_QUAL=20
FS_MASK=30

# MaskPrimers run parameters
MP_R1_MODE="mask"
MP_R2_MODE="cut"
MP_R1_MAXLEN=100
MP_R2_MAXLEN=100
MP_R1_MAXERR=0.2
MP_R2_MAXERR=0.2

# AssemblePairs-align run parameters
AP_ALN_SCANREV=true
AP_ALN_MAXERR=0.3
AP_ALN_MINLEN=8
AP_ALN_ALPHA=1e-5

# AssemblePairs-reference run parameters
AP_REF_MINIDENT=0.5
AP_REF_EVALUE=1e-5
AP_REF_MAXHITS=100
REF_FILE="/scratch2/kleinstein/germlines/IMGT_Human_IGV_ungapped_2014-08-23.fasta"
#REF_FILE="/scratch2/kleinstein/germlines/IMGT_Mouse_IGV_ungapped_2014-11-22.fasta"
USEARCH_EXEC=$HOME/bin/usearch

# CollapseSeq run parameters
CS_KEEP=true
CS_MISS=0

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

# Assemble paired ends via mate-pair alignment
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"

if $AP_ALN_SCANREV; then
    AssemblePairs.py align -1 $R1_FILE -2 $R2_FILE --coord illumina --rc tail \
        --minlen $AP_ALN_MINLEN --maxerror $AP_ALN_MAXERR --alpha $AP_ALN_ALPHA \
        --scanrev --failed --outname "${OUTNAME}-ALN" --outdir . \
        --log AssembleAlignLog.log --nproc $NPROC >> $PIPELINE_LOG 2> $ERROR_LOG
else
    AssemblePairs.py align -1 $R1_FILE -2 $R2_FILE --coord illumina --rc tail \
        --minlen $AP_ALN_MINLEN --maxerror $AP_ALN_MAXERR --alpha $AP_ALN_ALPHA \
        --failed --outname "${OUTNAME}-ALN" --outdir . --log AssembleAlignLog.log \
        --nproc $NPROC >> $PIPELINE_LOG 2> $ERROR_LOG
fi

# Assemble paired ends via alignment against V-region reference database
if $REFERENCE_ASSEMBLY; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs reference"
    AssemblePairs.py reference -1 "${OUTNAME}-ALN-1_assemble-fail.fastq" \
        -2 "${OUTNAME}-ALN-2_assemble-fail.fastq" -r $REF_FILE --coord illumina \
        --minident $AP_REF_MINIDENT --evalue $AP_REF_EVALUE --maxhits $AP_REF_MAXHITS \
        --outname "${OUTNAME}-REF" --log AssembleRefLog.log --exec $USEARCH_EXEC \
        --nproc $NPROC >> $PIPELINE_LOG 2> $ERROR_LOG
    cat "${OUTNAME}-ALN_assemble-pass.fastq" "${OUTNAME}-REF_assemble-pass.fastq" > \
        "${OUTNAME}-CAT_assemble-pass.fastq"
    FS_FILE="${OUTNAME}-CAT_assemble-pass.fastq"
else
    FS_FILE="${OUTNAME}-ALN_assemble-pass.fastq"
fi

# Remove low quality reads
if $FILTER_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
    FilterSeq.py quality -s $FS_FILE -q $FS_QUAL --outname "${OUTNAME}" \
        --log QualityLog.log --nproc $NPROC >> $PIPELINE_LOG  2> $ERROR_LOG
    MP_FILE="${OUTNAME}_quality-pass.fastq"
else
    MP_FILE=$FS_FILE
fi

# Identify and mask V-region primers
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align (R1)"
MaskPrimers.py align -s $MP_FILE -p $R1_PRIMERS --mode $MP_R1_MODE \
    --maxlen $MP_R1_MAXLEN --maxerror $MP_R1_MAXERR --outname "${OUTNAME}-PR1" \
    --log PrimerLogR1.log --nproc $NPROC >> $PIPELINE_LOG 2> $ERROR_LOG

# Rename V-region primer field
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders rename (R1)"
ParseHeaders.py rename -s "${OUTNAME}-PR1_primers-pass.fastq" \
    -f PRIMER -k VPRIMER --outname "${OUTNAME}-PR1" > /dev/null

# Identify and mask C-region primers
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align (R2)"
MaskPrimers.py align -s "${OUTNAME}-PR1_reheader.fastq" -p $R2_PRIMERS --mode $MP_R2_MODE \
    --maxlen $MP_R2_MAXLEN --maxerror $MP_R2_MAXERR --revpr --skiprc --outname "${OUTNAME}-PR2" \
    --log PrimerLogR2.log --nproc $NPROC >> $PIPELINE_LOG 2> $ERROR_LOG

# Rename C-region primer field
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders rename (R2)"
ParseHeaders.py rename -s "${OUTNAME}-PR2_primers-pass.fastq" \
    -f PRIMER -k CPRIMER --outname "${OUTNAME}-PR2" > /dev/null

# Mask low quality positions
if $MASK_LOWQUAL; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq maskqual"
    FilterSeq.py maskqual -s "${OUTNAME}-PR2_reheader.fastq" -q $FS_MASK \
        --outname "${OUTNAME}-FIN" --log MaskqualLog.log \
        --nproc $NPROC >> $PIPELINE_LOG 2> $ERROR_LOG
    CS_FILE="${OUTNAME}-FIN_maskqual-pass.fastq"
else
    CS_FILE="${OUTNAME}-PR2_reheader.fastq"
fi

# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
if $CS_KEEP; then
    CollapseSeq.py -s $CS_FILE -n $CS_MISS --uf CPRIMER --cf VPRIMER --act set \
        --inner --keepmiss --outname "${OUTNAME}-FIN" >> $PIPELINE_LOG 2> $ERROR_LOG
else
    CollapseSeq.py -s $CS_FILE -n $CS_MISS --uf CPRIMER --cf VPRIMER --act set \
    --inner --outname "${OUTNAME}-FIN" >> $PIPELINE_LOG 2> $ERROR_LOG
fi

# Filter to sequences with at least 2 supporting sources
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s "${OUTNAME}-FIN_collapse-unique.fastq" -f DUPCOUNT --num 2 \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Create table of final repertoire
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
ParseHeaders.py table -s "${OUTNAME}-FIN_collapse-unique.fastq" \
    -f ID CPRIMER VPRIMER DUPCOUNT --outname "Unique" >> $PIPELINE_LOG 2> $ERROR_LOG
ParseHeaders.py table -s "${OUTNAME}-FIN_collapse-unique_atleast-2.fastq" \
    -f ID CPRIMER VPRIMER DUPCOUNT --outname "Final-Unique-Atleast2" \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l AssembleAlignLog.log -f ID LENGTH OVERLAP ERROR PVALUE FIELDS1 FIELDS2 \
    > /dev/null  2> $ERROR_LOG &
ParseLog.py -l PrimerLogR[1-2].log -f ID BARCODE PRIMER ERROR \
    > /dev/null  2> $ERROR_LOG &
if $REFERENCE_ASSEMBLY; then
    ParseLog.py -l AssembleRefLog.log -f ID REFID LENGTH OVERLAP GAP EVALUE1 EVALUE2 IDENTITY FIELDS1 FIELDS2 \
    > /dev/null  2> $ERROR_LOG &
fi
if $FILTER_LOWQUAL; then
    ParseLog.py -l QualityLog.log -f ID QUALITY > /dev/null &
fi
if $MASK_LOWQUAL; then
    ParseLog.py -l MaskqualLog.log -f ID MASKED > /dev/null  2> $ERROR_LOG &
fi
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    LOG_FILES_ZIP=$(ls *LogR[1-2].log *Log.log)
    tar -zcf LogFiles.tar $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP

    TEMP_FILES_ZIP=$(ls *.fastq | grep -vP "collapse-unique.fastq|atleast-2.fastq")
    tar -zcf TempFiles.tar $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
fi

# End
printf "DONE\n\n"
cd ../

