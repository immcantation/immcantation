#!/usr/bin/env bash
# Super script to run a simplified pRESTO 0.5.1, IgBLAST and Change-O 0.3.1 pipeline
#
# Author:  Jason Anthony Vander Heiden
# Date:    2016.02.03
#
# Required Arguments:
#   $1 = FASTA sequence file
#   $2 = V segment primer file
#   $3 = J segment primer file
#   $4 = output directory
#   $5 = output file prefix
#   $6 = number of subprocesses for multiprocessing tools

# Capture command line parameters
INPUT_FILE=$(readlink -f $1)
V_PRIMERS=$(readlink -f $2)
J_PRIMERS=$(readlink -f $3)
OUTDIR=$4
OUTNAME=$5
NPROC=$6

# Define pipeline steps
ZIP_FILES=true
PARSE_LOGS=true

# Universal parameters
REF_GAPPED="${HOME}/workspace/igpipeline/germlines/IMGT/Human/VDJ"

# MaskPrimers run parameters
MP_V_MODE="mask"
MP_J_MODE="mask"
MP_V_START=0
MP_J_START=0
MP_V_MAXERR=0.2
MP_J_MAXERR=0.2
MP_V_MAXLEN=30
MP_J_MAXLEN=30

# CollapseSeq run parameters
CS_KEEP=false
CS_MISS=20

# DefineClones run parameters
DC_MODEL=hs1f
DC_DIST=0.2
DC_ACT=first

# Create germlines parameters
CG_GERM=dmask
CG_SFIELD=SEQUENCE_IMGT
CG_VFIELD=V_CALL

# Define IgBLAST command
IGDATA="${HOME}/apps/igblast-1.4.0"
IGBLAST_DB="${IGDATA}/database"
IGBLAST_CMD="${IGDATA}/igblastn \
    -germline_db_V ${IGBLAST_DB}/imgt_human_IG_V \
    -germline_db_D ${IGBLAST_DB}/imgt_human_IG_D \
    -germline_db_J ${IGBLAST_DB}/imgt_human_IG_J \
    -auxiliary_data ${IGDATA}/optional_file/human_gl.aux \
    -domain_system imgt -ig_seqtype Ig -organism human \
    -outfmt '7 std qseq sseq btop'"

# Define log files
PIPELINE_LOG="Pipeline.log"
ERROR_LOG="Pipeline.err"

# Make output directory and empty log files
mkdir -p $OUTDIR; cd $OUTDIR
echo '' > $PIPELINE_LOG
echo '' > $ERROR_LOG

# Start
echo "DIRECTORY: ${OUTDIR}"
echo "PRESTO VERSIONS:"
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
echo "CHANGEO VERSIONS:"
echo "  $(CreateGermlines.py --version 2>&1)"
echo "  $(DefineClones.py --version 2>&1)"
echo "  $(MakeDb.py --version 2>&1)"
echo "  $(ParseDb.py --version 2>&1)"
echo "IGBLAST VERSION:"
echo "  $(${IGBLAST_CMD} -version | grep 'Package' |sed s/'Package: '//)"
echo -e "\nSTART"
STEP=0

# Identify primers
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s $INPUT_FILE -p $V_PRIMERS --mode $MP_V_MODE \
    --maxlen $MP_V_MAXLEN --maxerror $MP_V_MAXERR \
    --nproc $NPROC --log PrimerVLog.log --outname "${OUTNAME}-V" --outdir . \
    >> $PIPELINE_LOG 2> $ERROR_LOG
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s "${OUTNAME}-V_primers-pass.fasta" -p $J_PRIMERS --mode $MP_J_MODE \
    --maxlen $MP_J_MAXLEN --maxerror $MP_J_MAXERR --revpr --skiprc \
    --nproc $NPROC --log PrimerJLog.log --outname "${OUTNAME}-J" \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Expand primer field
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders expand"
ParseHeaders.py expand -s "${OUTNAME}-J_primers-pass.fasta" \
    -f PRIMER --outname "${OUTNAME}-RH1" \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Rename primer fields
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders rename"
ParseHeaders.py rename -s "${OUTNAME}-RH1_reheader.fasta" \
    -f PRIMER1 PRIMER2 -k VPRIMER CPRIMER --outname "${OUTNAME}-RH2" \
    >> $PIPELINE_LOG 2> $ERROR_LOG

# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
if $CS_KEEP; then
    CollapseSeq.py -s "${OUTNAME}-RH2_reheader.fasta" -n $CS_MISS --cf VPRIMER --act set \
    --inner --keepmiss --outname "${OUTNAME}" >> $PIPELINE_LOG 2> $ERROR_LOG
else
    CollapseSeq.py -s "${OUTNAME}-RH2_reheader.fasta" -n $CS_MISS --cf VPRIMER --act set \
    --inner --outname "${OUTNAME}" >> $PIPELINE_LOG 2> $ERROR_LOG
fi

# Align V(D)J segments using IgBLAST
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "igblastn"
IGBLAST_RUN="${IGBLAST_CMD} -query ${OUTNAME}_collapse-unique.fasta \
    -out ${OUTNAME}_collapse-unique.fmt7 \
    -num_threads ${NPROC}"
echo -e "   START> igblastn" >> $PIPELINE_LOG
echo -e "    FILE> ${OUTNAME}_collapse-unique.fasta\n" >> $PIPELINE_LOG
echo -e "PROGRESS> [Running]" >> $PIPELINE_LOG
eval $IGBLAST_RUN >> $PIPELINE_LOG 2> $ERROR_LOG
echo -e "PROGRESS> [Done   ]\n" >> $PIPELINE_LOG
echo -e "  OUTPUT> ${OUTNAME}_collapse-unique.fmt7" >> $PIPELINE_LOG
echo -e "     END> igblastn\n" >> $PIPELINE_LOG

# Parse IgBLAST output
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MakeDb igblast"
MakeDb.py igblast -i  "${OUTNAME}_collapse-unique.fmt7" \
    -s  "${OUTNAME}_collapse-unique.fasta" -r $REF_GAPPED \
    --scores --regions --outname "${OUTNAME}" \
    >> $PIPELINE_LOG 2> $ERROR_LOG

printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseDb select"
ParseDb.py select -d "${OUTNAME}_db-pass.tab" -f FUNCTIONAL -u T \
    --outname "${OUTNAME}" >> $PIPELINE_LOG 2> $ERROR_LOG

# Assign clones
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "DefineClones bygroup"
DefineClones.py bygroup -d "${OUTNAME}_parse-select.tab" --model $DC_MODEL \
    --dist $DC_DIST --mode gene --act $DC_ACT --nproc $NPROC --outname "${OUTNAME}" \
    --log CloneLog.log >> $PIPELINE_LOG 2> $ERROR_LOG

# Create germlines
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CreateGermlines"
CreateGermlines.py -d "${OUTNAME}_clone-pass.tab" -r $REF_GAPPED -g $CG_GERM \
    --sf $CG_SFIELD --vf $CG_VFIELD --cloned --outname "${OUTNAME}" \
    --log GermLog.log >> $PIPELINE_LOG 2> $ERROR_LOG

# Process log files
if $PARSE_LOGS; then
    # Create table of final repertoire
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
    ParseHeaders.py table -s "${OUTNAME}_collapse-unique.fasta" \
        -f ID VPRIMER JPRIMER DUPCOUNT --outname "Unique" \
        >> $PIPELINE_LOG 2> $ERROR_LOG

    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
    ParseLog.py -l Primer[VJ]Log.log -f ID PRIMER ERROR \
        > /dev/null  2> $ERROR_LOG &
    ParseLog.py -l CloneLog.log -f VALLELE DALLELE JALLELE JUNCLEN SEQUENCES CLONES \
        > /dev/null  2> $ERROR_LOG &
fi

# Zip intermediate and log files
if $ZIP_FILES; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Zipping files"
    LOG_FILES_ZIP=$(ls *Log.log)
    tar -cf LogFiles.tar $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP
    gzip LogFiles.tar

    TEMP_FILES_ZIP=$(ls *.tab *.fasta | grep -v "collapse-unique\|germ-pass.tab\|table.tab\|headers.tab")
    tar -cf TempFiles.tar $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
    gzip TempFiles.tar
fi

# End
echo -e "DONE\n"
cd ..