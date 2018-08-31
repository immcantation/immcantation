#!/usr/bin/env bats

# Run parametes
DATE=$(date +"%Y.%m.%d")
DATA_DIR="data"
DATA_DIR=$(readlink -f ${DATA_DIR})
RUN_DIR="run/changeo-${DATE}"

# Run parameters
NPROC=2
OUTDIR=false
GERMLINES="${HOME}/share/imgt/human/vdj"
V_GERMLINES="${HOME}/share/igblast/fasta/imgt_human_ig_v.fasta"
#GERMLINES="/usr/local/share/imgt/human/vdj"
#FORMAT=""
#FORMAT="--format airr"
#FASTA="--fasta"
#CLUSTER="cd-hit-est"
ALIGNER="blastn"
#FAILED=""
FAILED="--failed"

# Create output parent
mkdir -p ${RUN_DIR}/logs ${RUN_DIR}/console ${RUN_DIR}/output
RUN_DIR=$(readlink -f ${RUN_DIR})

# Create output parent
mkdir -p ${RUN_DIR}

# Changeo
AlignRecords across -d S43_IGBLAST_AIRR_db-pass.tsv --sf sequence junction_nt --calls v d j --format airr $OUTPUT $LOG $FAILED $NPROC
AlignRecords across -d RQ2860-2_MYG91N_db-pass_germ-pass.tab --sf SEQUENCE_INPUT GERMLINE_IMGT_D_MASK JUNCTION --gf JUNCTION_LENGTH \
    --calls v d j --nproc 2 $OUTPUT $LOG $FAILED
AlignRecords block -d S43_IGBLAST_AIRR_db-pass.tsv --sf sequence junction_nt --gf junction_length --calls v d j --format airr \
    $OUTPUT $LOG $FAILED $NPROC
AlignRecords block -d RQ2410_SR_Norm_CHAIN=heavy_clone-pass_germ-pass.tab --sf SEQUENCE GERMLINE_GAP_D_MASK JUNCTION --gf CLONE --calls v d j \
    $OUTPUT $LOG $FAILED $NPROC
AlignRecords within -d S43_IGBLAST_AIRR_db-pass.tsv --sf sequence junction --format airr $OUTPUT $LOG $FAILED $NPROC
AlignRecords within -d RQ2860-2_MYG91N_db-pass_germ-pass.tab --sf SEQUENCE_INPUT GERMLINE_IMGT_D_MASK JUNCTION $OUTPUT $LOG $FAILED $NPROC
AssignGenes igblast -s S43_atleast-2.fasta -b $USER_HOME$/share/igblast --format blast --organism human $OUTPUT 
ConvertDb airr -d SUBJECT-AR02_SORT-Memory.tsv $OUTPUT
ConvertDb baseline -d SUBJECT-AR02_SORT-Memory.tsv --sf SEQUENCE_IMGT --gf GERMLINE_IMGT_D_MASK --mf V_CALL J_CALL --cf CLONE $OUTPUT
ConvertDb changeo -d SRR1383447_airr_N5K.tsv $OUTPUT
ConvertDb fasta -d S43_IGBLAST_CHANGEO_db-pass.tsv --sf SEQUENCE_IMGT --if SEQUENCE_ID --mf V_CALL J_CALL $OUTPUT
ConvertDb genbank -d S43_IGBLAST_CHANGEO_db-pass.tsv --inf IgBLAST:1.7.0 --organism "Homo sapiens" --sex Male --tissue "Peripheral blood" \
    --cf CPRIMER --nf DUPCOUNT --if MID --asis-id -y genbank_test.yaml --asn --sbt output/template.sbt $OUTPUT
ConvertDb genbank -d SUBJECT-AR02_SORT-Memory.tsv --inf IMGT/HighV-QUEST:1.5.5 --organism "Homo sapiens" --cf CREGION \
    --cell-type "memory B cell" --isolate AR03 --asis-id $OUTPUT
CreateGermlines -d S43_IGBLAST_AIRR_db-pass_clone-pass.tsv -r germlines_gapped -g vonly dmask full regions --format airr $LOG $OUTPUT
CreateGermlines -d S43_IGBLAST_AIRR_db-pass_clone-pass.tsv -r germlines_gapped -g vonly dmask full regions --cloned --format airr $LOG $OUTPUT
CreateGermlines -d AA_RMHBx1_db-pass_clone-pass.tab --cloned -r $USER_HOME$/share/imgt/human/vdj -g vonly dmask full regions $FAILED $LOG $OUTPUT
CreateGermlines -d S43_IGBLAST_CHANGEO_db-pass.tab -r $USER_HOME$/share/imgt/human/vdj -g vonly dmask full regions $LOG $FAILED $OUTPUT
DefineClones -d S43_IGBLAST_AIRR_db-pass.tsv --model ham --dist 0.15 --mode gene --maxmiss 5 --act first --format airr $FAILED $LOG $OUTPUT
DefineClones -d A79HP_MG41M_db-pass.tab --model ham --dist 0.15 --mode gene --maxmiss 0 --act set $FAILED $OUTPUT $LOG
MakeDb igblast -i igblast/v1.7/S43_atleast-2.fmt7 -s igblast/v1.7/S43_atleast-2.fasta -r $PROJECT_DIR$/../germlines/IMGT/human/vdj \
    --scores --partial --regions --format airr $OUTPUT
MakeDb igblast -i igblast/v1.7/S43_atleast-2.fmt7 -s igblast/v1.7/S43_atleast-2.fasta -r $USER_HOME$/share/imgt/human/vdj \
    --regions --scores --cdr3 $OUTPUT $FAILED
MakeDb ihmm -i ihmmune-align/Short_ihmm.txt -s ihmmune-align/ShortSimSeqs_sequences.fasta -r $PROJECT_DIR$/../germlines/IMGT/human/vdj \
    --scores --regions --partial $OUTPUT $FAILED
MakeDb imgt -i imgt/S43_atleast-2.txz -s imgt/S43_atleast-2.fasta --regions --scores --junction --format airr $FAILED $OUTPUT
MakeDb imgt -i imgt/HD23_Naive.txz -s imgt/HD23_Naive-FIN_collapse-unique_atleast-2.fasta \
    -r $PROJECT_DIR$/../germlines/IMGT/human/vdj/imgt_human_IGHV.fasta $PROJECT_DIR$/../germlines/IMGT/human/vdj/imgt_human_IGHD.fasta $PROJECT_DIR$/../germlines/IMGT/human/vdj/imgt_human_IGHJ.fasta \
    --regions --scores --junction $FAILED $OUTPUT
ParseDb add -d GL10870_nodes.tab -f ADD_1 ADD_2 -u 1 2 $OUTPUT
ParseDb delete -d GL10870_nodes.tab -f V_CALL -u "Inferred|Germline" --regex $OUTPUT
ParseDb drop -d GL10870_nodes.tab -f V_CALL J_CALL $OUTPUT
ParseDb index -d GL10870_nodes.tab -f INDEX $OUTPUT
ParseDb rename -d GL10870_nodes.tab -f V_CALL J_CALL -k V_CALL_RENAME J_CALL_RENAME $OUTPUT
ParseDb select -d GL10870_nodes.tab -f V_CALL J_CALL -u IGH --logic all --regex $OUTPUT
ParseDb sort -d S43_IMGT_CHANGEO_db-pass.tsv -f DUPCOUNT --num --descend $OUTPUT
ParseDb split -d S43_IMGT_CHANGEO_db-pass.tsv -f CPRIMER $OUTPUT
ParseDb update -d GL10870_nodes.tab -f SAMPLE_TYPE -u Inferred Blood -t Test Test $OUTPUT
