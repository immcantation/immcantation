#!/usr/bin/env bash

# Run parameters
# sudo singularity build immcantation-devel.img docker-daemon://immcantation/suite:devel
# apptainer build immcantation-devel.img docker-daemon://immcantation/suite:devel
VERSION="devel"
IMAGE=$(realpath immcantation-devel.img)
DATE=$(date +"%Y.%m.%d")
DATA_DIR=$(realpath data)
RUN_DIR="run/singularity-${DATE}-${VERSION}"
SAMPLE=HD13M
NPROC=8
EXT="tsv"

# Create output parent
mkdir -p ${RUN_DIR}
RUN_DIR=$(realpath ${RUN_DIR})

# PhiX
READS_R1=/data/sequences/HD13M_10K_R1.fastq
OUT_DIR="/scratch/phix"

singularity exec -e -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE preprocess-phix \
    -s $READS_R1 -o $OUT_DIR -p $NPROC

# AbSeq
READS_R1=/data/sequences/HD13M_10K_R1.fastq
READS_R2=/data/sequences/HD13M_10K_R2.fastq
YAML=/data/report.yaml
OUT_DIR="/scratch/presto"

singularity exec -e -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE presto-abseq \
    -1 $READS_R1 -2 $READS_R2 -y $YAML -n $SAMPLE -o $OUT_DIR -p $NPROC

# Clontech
SAMPLE=HD13M
READS_R1=/data/sequences/HD13M_10K_R1.fastq
READS_R2=/data/sequences/HD13M_10K_R2.fastq
CREGION=/usr/local/share/protocols/Universal/Human_IG_CRegion_RC.fasta
VREF=/usr/local/share/igblast/fasta/imgt_human_ig_v.fasta
YAML=/data/report.yaml
OUT_DIR="/scratch/clontech"

singularity exec -e -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE \
        presto-clontech -1 $READS_R1 -2 $READS_R2 -j $CREGION -r $VREF \
        -y $YAML -n $SAMPLE -o $OUT_DIR -p $NPROC


# Clontech-umi
SAMPLE=HD13M
READS_R1=/data/sequences/HD13M_10K_R1.fastq
READS_R2=/data/sequences/HD13M_10K_R2.fastq
CREGION=/usr/local/share/protocols/Universal/Human_IG_CRegion_RC.fasta
VREF=/usr/local/share/igblast/fasta/imgt_human_ig_v.fasta
OUT_DIR="/scratch/clontech-umi"

singularity exec -e -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE \
        presto-clontech-umi -1 $READS_R1 -2 $READS_R2 -j $CREGION -r $VREF \
        -n $SAMPLE -o $OUT_DIR -p $NPROC

# 10X
SAMPLE=PBMC2B
READS=/data/sequences/PBMC2B.fasta
ANNOTATIONS=/data/sequences/PBMC2B_annotations.csv
DIST=0.15
OUT_DIR="/scratch/10x"

singularity exec -e -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE \
        changeo-10x -s $READS -a $ANNOTATIONS -x $DIST \
        -n $SAMPLE -o $OUT_DIR -p $NPROC -z


# IgBLAST
SAMPLE=HD13M
READS="/scratch/presto/${SAMPLE}-final_collapse-unique_atleast-2.fastq"
OUT_DIR="/scratch/changeo"

singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE changeo-igblast \
    -s $READS -n $SAMPLE -o $OUT_DIR -p $NPROC

# Change-O cloning
SAMPLE=HD13M
DB="/scratch/changeo/${SAMPLE}_db-pass.${EXT}"
OUT_DIR="/scratch/changeo"
DIST=0.15

singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE  \
    changeo-clone -d $DB -x $DIST -n $SAMPLE -o $OUT_DIR -p $NPROC -z


# TIgGER
SAMPLE=HD13M
DB="/scratch/changeo/${SAMPLE}_db-pass.${EXT}"
V_FIELD="v_call_genotyped"
MINSEQ=2
MINGERM=2
FIND_UNMUTATED=F
OUT_DIR="/scratch/changeo"
DB_H="${SAMPLE}_igh_db-pass"

# Subset to heavy chain
singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE \
        ParseDb.py select -d $DB -f v_call -u IGH --regex \
        --outname "${DB_H}" --outdir $OUT_DIR

singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE tigger-genotype \
        -d $OUT_DIR/$DB_H'_parse-select.tsv' -v $V_FIELD \
        -x $MINSEQ -y $MINGERM -u $FIND_UNMUTATED \
        -n $SAMPLE -o $OUT_DIR -p $NPROC

# SHazaM threshold
SAMPLE=HD13M
DB="/scratch/changeo/${SAMPLE}_db-pass.${EXT}"
OUT_DIR="/scratch/changeo"

singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE \
        shazam-threshold -d $DB -m density \
        -n "${SAMPLE}-1" -o $OUT_DIR -p $NPROC

singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE \
        shazam-threshold -d $DB -m gmm \
        -n "${SAMPLE}-2" -o $OUT_DIR -p $NPROC

singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE \
        shazam-threshold -d $DB --subsample 100 --repeats 2 \
        -n "${SAMPLE}-3" -o $OUT_DIR -p $NPROC

singularity exec -B $DATA_DIR:/data -B $RUN_DIR:/scratch $IMAGE \
          shazam-threshold -d $DB -m none \
          -n "${SAMPLE}-4" -o $OUT_DIR -p $NPROC
