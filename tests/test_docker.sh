#!/usr/bin/env bats

# Run parameters
DATE=$(date +"%Y.%m.%d")
IMAGE=kleinstein/immcantation:devel
DATA_DIR="/home/jason/workspace/igpipeline/immcantation/tests/data"
RUN_DIR="/home/jason/workspace/igpipeline/immcantation/tests/run/${DATE}"
SAMPLE=HD13M
NPROC=2
EXT="tab"

# Create output parent
mkdir -p $RUN_DIR

# PhiX
@test "preprocess-phix" {
	READS_R1=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq
	READS_R2=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R2_001.fastq
	OUT_DIR="/scratch/phix/${SAMPLE}"
	LOG="/scratch/run_preprocess-phix.out"

	docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		preprocess-phix -s $READS_R1 -o $OUT_DIR -p $NPROC \
		2>&1 $LOG
}
	
# pRESTO
@test "presto-abseq" {
	READS_R1=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq
	READS_R2=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R2_001.fastq
	YAML=/data/report.yaml
	OUT_DIR="/scratch/presto/${SAMPLE}"
	LOG="/scratch/run_presto-abseq.out"

	docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		presto-abseq -1 $READS_R1 -2 $READS_R2 -y $YAML -n $SAMPLE -o $OUT_DIR -p $NPROC \
		2>&1 $LOG
}

# IgBLAST
@test "changeo-igblast" {
	READS="/scratch/presto/${SAMPLE}/${SAMPLE}-final_collapse-unique_atleast-2.fastq"
	OUT_DIR="/scratch/changeo/${SAMPLE}"
	LOG="/scratch/run_changeo-igblast.out"

	docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		changeo-igblast -s $READS -n $SAMPLE -o $OUT_DIR -p $NPROC \
		2>&1 $LOG
}

# TIgGER
@test "tigger-genotype" {
	DB="/scratch/changeo/${SAMPLE}/${SAMPLE}_db-pass.${EXT}"
	OUT_DIR="/scratch/changeo/${SAMPLE}"
	LOG="/scratch/run_tigger-genotype.out"

	docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		tigger-genotype -d $DB -n $SAMPLE -o $OUT_DIR -p $NPROC \
		2>&1 $LOG
}

# SHazaM threshold
@test "shazam-threshold" {
	DB="/scratch/changeo/${SAMPLE}/${SAMPLE}_genotyped.${EXT}"
	OUT_DIR="/scratch/changeo/${SAMPLE}"
	LOG="/scratch/run_shazam-threshold.out"

	docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		shazam-threshold -d $DB -n $SAMPLE -o $OUT_DIR -p $NPROC \
		2>&1 $LOG
}

# Change-O cloning
@test "changeo-clone" {
	DB="/scratch/changeo/${SAMPLE}/${SAMPLE}_genotyped.${EXT}"
	OUT_DIR="/scratch/changeo/${SAMPLE}"
	DIST=0.15
	LOG="/scratch/run_changeo-clone.out"

	docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		changeo-clone -d $DB -x $DIST -n $SAMPLE -o $OUT_DIR -p $NPROC \
		2>&1 $LOG
}
