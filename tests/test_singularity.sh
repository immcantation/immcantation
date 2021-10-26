#!/usr/bin/env bash

# Run parameters
IMAGE=$(realpath immcantation-devel.img)
DATE=$(date +"%Y.%m.%d")
DATA_DIR=$(realpath data)
RUN_DIR="run/${DATE}"
SAMPLE=HD13M
NPROC=8
EXT="tab"

# Create output parent
mkdir -p ${RUN_DIR}
RUN_DIR=$(realpath ${RUN_DIR})

# PhiX
READS_R1=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq
READS_R2=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R2_001.fastq
OUT_DIR="/scratch/phix"

singularity exec -e -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE preprocess-phix \
    -s $READS_R1 -o $OUT_DIR -p $NPROC

# pRESTO
READS_R1=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq
READS_R2=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R2_001.fastq
YAML=/data/report.yaml
OUT_DIR="/scratch/presto"

singularity exec -e -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE presto-abseq \
    -1 $READS_R1 -2 $READS_R2 -y $YAML -n $SAMPLE -o $OUT_DIR -p $NPROC

# IgBLAST
READS="/scratch/presto/${SAMPLE}-final_collapse-unique_atleast-2.fastq"
OUT_DIR="/scratch/changeo"

singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE changeo-igblast \
    -s $READS -n $SAMPLE -o $OUT_DIR -p $NPROC

# TIgGER
DB="/scratch/changeo/${SAMPLE}_db-pass.${EXT}"
OUT_DIR="/scratch/changeo"

singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE tigger-genotype \
    -d $DB -n $SAMPLE -o $OUT_DIR -p $NPROC

# SHazaM threshold
DB="/scratch/changeo/${SAMPLE}_genotyped.${EXT}"
OUT_DIR="/scratch/changeo"

singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE shazam-threshold \
    -d $DB -n $SAMPLE -o $OUT_DIR -p $NPROC

# Change-O clones
DB="/scratch/changeo/${SAMPLE}_genotyped.${EXT}"
OUT_DIR="/scratch/changeo"
DIST=0.15

singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE changeo-clone \
    -d $DB -x $DIST -n $SAMPLE -o $OUT_DIR -p $NPROC
