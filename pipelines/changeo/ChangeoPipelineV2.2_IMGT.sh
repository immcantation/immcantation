#!/bin/bash
# Super script to run the Change-O clonal assignment pipeline on IMGT data 
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2015.04.07
# 
# Required Arguments:
#   $1 = IMGT zip file 
#   $2 = FASTA file that was submitted to IMGT
#   $3 = the folder containing the IMGT reference database sequences
#   $4 = output directory
#   $5 = output file prefix
#   $6 = number of subprocesses for multiprocessing tools


# Capture command line parameters
IMGT_FILE=$(readlink -f $1)
SEQ_FILE=$(readlink -f $2)
GERM_DIR=$(readlink -f $3)
OUTDIR=$(readlink -f $4)
OUTNAME=$5
NPROC=$6

# Define run parameters
LOG_RUNTIMES=true
ZIP_FILES=true
DEFINE_CLONES=true

# DefineClones parameters
DC_MODEL=m1n
DC_DIST=5
DC_ACT=first

# Create germlines parameters
CG_GERM=dmask
CG_SFIELD=SEQUENCE_IMGT
CG_VFIELD=V_CALL

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
echo "  $(CreateGermlines.py --version 2>&1)"
echo "  $(DefineClones.py --version 2>&1)"
echo "  $(MakeDb.py --version 2>&1)"
echo "  $(ParseDb.py --version 2>&1)"
echo -e "\nSTART"
STEP=0

# Parse IMGT output
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MakeDb imgt"
#echo MakeDb.py imgt -z $IMGT_FILE -s $SEQ_FILE --outdir . --clean ">>" $PIPELINE_LOG
MakeDb.py imgt -i $IMGT_FILE -s $SEQ_FILE --outname "${OUTNAME}" \
    --outdir . >> $PIPELINE_LOG 2> $ERROR_LOG

printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseDb select"
ParseDb.py select -d "${OUTNAME}_db-pass.tab" -f FUNCTIONAL -u T \
    --outname "${OUTNAME}" >> $PIPELINE_LOG 2> $ERROR_LOG
ParseDb.py select -d "${OUTNAME}_select-pass.tab" -f V_CALL J_CALL -u IGH \
    --logic all --regex --outname "${OUTNAME}_heavy" >> $PIPELINE_LOG 2> $ERROR_LOG
#ParseDb.py select -d "${OUTNAME}_select-pass.tab" -f V_CALL -u "IG[LK]" --regex \
#    --outname "${OUTNAME}_light" >> $PIPELINE_LOG 2> $ERROR_LOG

# Assign clones
if $DEFINE_CLONES; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "DefineClones bygroup"
    DefineClones.py bygroup -d "${OUTNAME}_parse-select.tab" --model $DC_MODEL \
        --dist $DC_DIST --mode gene --act $DC_ACT --nproc $NPROC --outname "${OUTNAME}" \
        --log CloneLog.log >> $PIPELINE_LOG 2> $ERROR_LOG
    CG_FILE="${OUTNAME}_clone-pass.tab"
else
    CG_FILE="${OUTNAME}_parse-select.tab"
fi

# Create germlines
if $DEFINE_CLONES; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CreateGermlines"
    CreateGermlines.py -d $CG_FILE -r $GERM_DIR -g $CG_GERM --sf $CG_SFIELD \
    --vf $CG_VFIELD --cloned --outname "${OUTNAME}" \
    --log GermLog.log >> $PIPELINE_LOG 2> $ERROR_LOG
else
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CreateGermlines"
    CreateGermlines.py -d $CG_FILE -r $GERM_DIR -g $CG_GERM --sf $CG_SFIELD \
        --vf $CG_VFIELD --outname "${OUTNAME}" \
        --log GermLog.log >> $PIPELINE_LOG 2> $ERROR_LOG
fi

# Zip intermediate and log files
if $ZIP_FILES; then
    LOG_FILES_ZIP=$(ls *Log.log)
    tar -cf LogFiles.tar $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP
    gzip LogFiles.tar

    TEMP_FILES_ZIP=$(ls *.tab | grep -v "db-pass.tab\|germ-pass.tab")
    tar -cf TempFiles.tar $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
    gzip TempFiles.tar
fi

# End
echo -e "DONE\n" 
cd ..