#!/bin/bash
# Run Change-O clonal assignment and germline creation
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2015.07.07
# 
# Required Arguments:
#   $1 = Change-O db file
#   $2 = the folder containing the IMGT reference database sequences
#   $3 = output directory
#   $4 = output filename prefix
#   $5 = number of subprocesses for multiprocessing tools


# Capture command line parameters
DB_FILE=$(readlink -f $1)
GERM_DIR=$(readlink -f $2)
OUTDIR=$3
OUTNAME=$4
NPROC=$5

# Define run parameters
ZIP_FILES=true

# DefineClones parameters
DC_MODEL=hs1f
DC_DIST=0.15
DC_MODE=gene
DC_ACT=set
DC_NORM=len
DC_LINK=single

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

# Assign clones
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "DefineClones bygroup"
DefineClones.py bygroup -d $DB_FILE --model $DC_MODEL --dist $DC_DIST \
    --mode $DC_MODE --act $DC_ACT --norm $DC_NORM --link $DC_LINK \
    --nproc $NPROC --failed --outname "${OUTNAME}" --outdir . \
    --log CloneLog.log >> $PIPELINE_LOG 2> $ERROR_LOG

# Create germlines
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CreateGermlines"
CreateGermlines.py -d "${OUTNAME}_clone-pass.tab" -r $GERM_DIR -g $CG_GERM \
    --sf $CG_SFIELD --vf $CG_VFIELD --cloned --failed \
    --log GermLog.log >> $PIPELINE_LOG 2> $ERROR_LOG

# Zip intermediate and log files
if $ZIP_FILES; then
    LOG_FILES_ZIP=$(ls *Log.log)
    tar -zcf LogFiles.tar.gz $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP

    TEMP_FILES_ZIP=$(ls *.tab | grep -v "germ-pass.tab")
    tar -zcf TempFiles.tar.gz $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
fi

# End
echo -e "DONE\n" 
cd ..