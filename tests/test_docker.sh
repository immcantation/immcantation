# Run parameters
DATE=$(date +"%Y.%m.%d")
IMAGE=kleinstein/immcantation:devel
DATA_DIR=/home/jason/workspace/igpipeline/immcantation/tests/data
SAMPLE_NAME=HD13M
NPROC=2

# Output parent
mkdir -p run/${DATE}

# PhiX
READS_R1=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq
READS_R2=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R2_001.fastq

docker run -v $DATA_DIR:/data:z $IMAGE preprocess-phix \
	-s $READS_R1 -o $OUT_DIR -p $NPROC \
	| tee run/${DATE}/run_phix_1.out
docker run -v $DATA_DIR:/data:z $IMAGE preprocess-phix \
	-s $READS_R2 -o $OUT_DIR -p $NPROC \
	| tee run/${DATE}/run_phix_2.out
	
# pRESTO
READS_R1=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq
READS_R2=/data/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R2_001.fastq
YAML=/data/report.yaml
OUT_DIR="/data/run/${DATE}/presto/${SAMPLE_NAME}"

docker run -v $DATA_DIR:/data:z $IMAGE presto-abseq \
    -1 $READS_R1 -2 $READS_R2 -y $YAML -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
    | tee run/${DATE}/run_presto.out

# IgBLAST
READS="/data/run/${DATE}/presto/${SAMPLE_NAME}/${SAMPLE_NAME}-final_collapse-unique_atleast-2.fastq"
OUT_DIR="/data/run/${DATE}/changeo/${SAMPLE_NAME}"

docker run -v $DATA_DIR:/data:z $IMAGE changeo-igblast \
    -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
    | tee run/${DATE}/run_igblast.out

# TIgGER
DB="/data/run/${DATE}/changeo/${SAMPLE_NAME}/${SAMPLE_NAME}_db-pass.tab"
OUT_DIR="/data/run/${DATE}/changeo/${SAMPLE_NAME}"

docker run -v $DATA_DIR:/data:z $IMAGE tigger-genotype \
    -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
    | tee run/${DATE}/run_genotype.out

# SHazaM threshold
DB="/data/run/${DATE}/changeo/${SAMPLE_NAME}/${SAMPLE_NAME}_genotyped.tab"
OUT_DIR="/data/run/${DATE}/changeo/${SAMPLE_NAME}"

docker run -v $DATA_DIR:/data:z $IMAGE shazam-threshold \
	-d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
	| tee run/${DATE}/run_threshold.out

# Change-O clones
DB="/data/run/${DATE}/changeo/${SAMPLE_NAME}/${SAMPLE_NAME}_genotyped.tab"
OUT_DIR="/data/run/${DATE}/changeo/${SAMPLE_NAME}"
DIST=0.15

docker run -v $DATA_DIR:/data:z $IMAGE changeo-clone \
	-d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
	| tee run/${DATE}/run_clone.out
