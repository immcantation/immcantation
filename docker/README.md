## Download

```
docker pull javh/immcantation:release
```

## List available pipelines

```
docker run -it javh/immcantation:release
```

## Specific pipline help

```
docker run -it javh/immcantation:release presto-abseq -h
```


## Invoke shell session within image

```
DATA_DIR=~/workspace/igpipeline/docker/tmp
docker run -it -v $DATA_DIR:/data:z javh/immcantation:release bash
```

## Run pRESTO pipeline for AbSeq v3 data

```
# Arguments
DATA_DIR=~/workspace/igpipeline/docker/tmp
READS_R1=/data/raw/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq
READS_R2=/data/raw/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R2_001.fastq
PRIMERS_R1=/data/primers/AbSeqV3_Human_R1CPrimers.fasta
PRIMERS_R2=/data/primers/AbSeqV3_Human_R2TSPrimers.fasta
CREGION=/data/primers/AbSeqV3_Human_InternalCRegion.fasta
YAML=/data/test.yaml
SAMPLE_NAME=AAYHL_HD13M
OUT_DIR=/data/presto/AAYHL_HD13M
NPROC=4

# Run pipeline in docker image.
# Note: mounting a host directory as a volume with the 'z' option is essential.
docker run -it -v $DATA_DIR:/data:z javh/immcantation:release presto-abseq \
    -1 $READS_R1 -2 $READS_R2 -j $PRIMERS_R1 -v $PRIMERS_R2 \
    -c $CREGION -y $YAML -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
    | tee run_presto.out
```

## Run IgBLAST pipeline

```
# Arguments
DATA_DIR=~/workspace/igpipeline/docker/tmp
READS=/data/presto/AAYHL_HD13M/AAYHL_HD13M-final_collapse-unique_atleast-2.fastq
SAMPLE_NAME=AAYHL_HD13M
OUT_DIR=/data/changeo/AAYHL_HD13M
NPROC=4

# Run pipeline in docker image.
# Note: mounting a host directory as a volume with the 'z' option is essential.
docker run -it -v $DATA_DIR:/data:z javh/immcantation:release changeo-igblast \
    -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
    | tee run_igblast.out
```

## Run clonal threshold inferrence

```
# Arguments
DATA_DIR=~/workspace/igpipeline/docker/tmp
DB=/data/changeo/AAYHL_HD13M/AAYHL_HD13M_db-pass.tab
SAMPLE_NAME=AAYHL_HD13M
OUT_DIR=/data/changeo/AAYHL_HD13M
NPROC=4

# Run pipeline in docker image.
# Note: mounting a host directory as a volume with the 'z' option is essential.
docker run -it -v $DATA_DIR:/data:z javh/immcantation:release shazam-threshold \
    -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
    | tee run_threshold.out
```

## Run clonal assignment pipeline

```
# Arguments
DATA_DIR=~/workspace/igpipeline/docker/tmp
DB=/data/changeo/AAYHL_HD13M/AAYHL_HD13M_db-pass.tab
DIST=0.15
SAMPLE_NAME=AAYHL_HD13M
OUT_DIR=/data/changeo/AAYHL_HD13M
NPROC=4

# Run pipeline in docker image.
# Note: mounting a host directory as a volume with the 'z' option is essential.
docker run -it -v $DATA_DIR:/data:z javh/immcantation:release changeo-clone \
    -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
    | tee run_clone.out
```

