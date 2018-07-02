#!/usr/bin/env bats

# Run parameters
VERSION="1.10.1"
IMAGE="kleinstein/immcantation:${VERSION}"
DATE=$(date +"%Y.%m.%d")
DATA_DIR=$(readlink -f data)
RUN_DIR="run/${DATE}"
SAMPLE=HD13M
NPROC=2
EXT="tab"

# Create output parent
mkdir -p ${RUN_DIR}
RUN_DIR=$(readlink -f ${RUN_DIR})

# PhiX
@test "preprocess-phix" {
	READS_R1=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq
	READS_R2=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R2_001.fastq
	OUT_DIR="/scratch/phix"

	run docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		preprocess-phix -s $READS_R1 -o $OUT_DIR -p $NPROC

	[ "$status" -eq 0 ]
}
	
# pRESTO
@test "presto-abseq" {
	READS_R1=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq
	READS_R2=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R2_001.fastq
	YAML=/data/report.yaml
	OUT_DIR="/scratch/presto"

	run docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		presto-abseq -1 $READS_R1 -2 $READS_R2 -y $YAML -n $SAMPLE -o $OUT_DIR -p $NPROC

	[ "$status" -eq 0 ]
}

# IgBLAST
@test "changeo-igblast" {
	READS="/scratch/presto/${SAMPLE}-final_collapse-unique_atleast-2.fastq"
	OUT_DIR="/scratch/changeo"

	run docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		changeo-igblast -s $READS -n $SAMPLE -o $OUT_DIR -p $NPROC

	[ "$status" -eq 0 ]
}

# TIgGER
@test "tigger-genotype" {
	DB="/scratch/changeo/${SAMPLE}_db-pass.${EXT}"
	OUT_DIR="/scratch/changeo"

	run docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		tigger-genotype -d $DB -n $SAMPLE -o $OUT_DIR -p $NPROC

	[ "$status" -eq 0 ]
}

# SHazaM threshold
@test "shazam-threshold" {
	DB="/scratch/changeo/${SAMPLE}_genotyped.${EXT}"
	OUT_DIR="/scratch/changeo"

	run docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		shazam-threshold -d $DB -n $SAMPLE -o $OUT_DIR -p $NPROC

	[ "$status" -eq 0 ]
}

# Change-O cloning
@test "changeo-clone" {
	DB="/scratch/changeo/${SAMPLE}_genotyped.${EXT}"
	OUT_DIR="/scratch/changeo"
	DIST=0.15

	run docker run -v $DATA_DIR:/data:z -v $RUN_DIR:/scratch:z $IMAGE \
		changeo-clone -d $DB -x $DIST -n $SAMPLE -o $OUT_DIR -p $NPROC

	[ "$status" -eq 0 ]
}
