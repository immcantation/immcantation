#!/usr/bin/env bash
# Super script to run IgBLAST and Change-O on 10x data
#
# Author:  Jason Anthony Vander Heiden, Ruoyi Jiang
# Date:    2019.05.15
#
# Arguments:
#   -s  FASTA or FASTQ sequence file.
#   -a  10x Genomics cellranger-vdj contig annotation file.
#       Must corresponding with the FASTA/FASTQ input file (all, filtered or consensus).
#   -r  Directory containing IMGT-gapped reference germlines.
#       Defaults to /usr/local/share/germlines/imgt/[species name]/vdj.
#   -g  Species name. One of human, mouse, rabbit, rat, or rhesus_monkey. Defaults to human.
#   -t  Receptor type. One of ig or tr. Defaults to ig.
#   -x  Distance threshold for clonal assignment. Specify 'auto' for automatic detection.
#       If unspecified, clonal assignment is not performed.
#   -m  Distance model for clonal assignment.
#       Defaults to the nucleotide Hamming distance model (nt). Options: nt, aa.
#   -e  Method to use for determining the optimal threshold. One of "gmm" or "density".
#       Defaults to "density".
#   -d  Curve fitting model. Applies only when method (-e) is 'gmm'. Options: are
#       'norm-norm', 'norm-gamma', 'gamma-norm' and 'gamma-gamma'.
#       Defaults to 'gamma-gamma'.
#   -u  Method to use for threshold selection. Applies only when method (-e) is 'gmm'.
#       One of 'optimal', 'intersect' and 'user'.
#       Defaults to 'user'.
#   -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file
#       and auxillary_data directories. Defaults to /usr/local/share/igblast.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the read 1 filename.
#   -o  Output directory. Will be created if it does not exist.
#       Defaults to a directory matching the sample identifier in the current working directory.
#   -f  Output format. One of changeo or airr. Defaults to airr.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -i  Specify to allow partial alignments.
#   -z  Specify to disable cleaning and compression of temporary files.
#   -h  Display help.

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -s  FASTA or FASTQ sequence file."
    echo -e "  -a  10x Genomics cellranger-vdj contig annotation CSV file.\n" \
            "     Must corresponding with the FASTA/FASTQ input file (all, filtered or consensus)."
    echo -e "  -r  Directory containing IMGT-gapped reference germlines.\n" \
            "     Defaults to /usr/local/share/germlines/imgt/[species name]/vdj."
    echo -e "  -g  Species name. One of human, mouse, rabbit, rat, or rhesus_monkey. Defaults to human."
    echo -e "  -t  Receptor type. One of ig or tr. Defaults to ig."
    echo -e "  -x  Distance threshold for clonal assignment. Specify \"auto\" for automatic detection.\n" \
            "     If unspecified, clonal assignment is not performed."
    echo -e "  -m  Distance model for clonal assignment.\n" \
            "     Defaults to the nucleotide Hamming distance model (ham)."
    echo -e "  -e Method to use for determining the optimal threshold. One of 'gmm' or 'density'. \n" \
            "     Defaults to 'density'."
    echo -e "  -d Curve fitting model. Applies only when method (-e) is 'gmm'. One of 'norm-norm',\n" \
            "      'norm-gamma', 'gamma-norm' and 'gamma-gamma'. \n" \
            "     Defaults to 'gamma-gamma'."
    echo -e "  -u Method to use for threshold selection. Applies only when method (-e) is 'gmm'. \n" \
            "     One of 'optimal', 'intersect' and 'user'. \n" \
            "     Defaults to 'user'."
    echo -e "  -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file\n" \
            "     and auxillary_data directories. Defaults to /usr/local/share/igblast."
    echo -e "  -n  Sample identifier which will be used as the output file prefix.\n" \
            "     Defaults to a truncated version of the sequence filename."
    echo -e "  -o  Output directory. Will be created if it does not exist.\n" \
            "     Defaults to a directory matching the sample identifier in the current working directory."
    echo -e "  -f  Output format. One of changeo or airr. Defaults to airr."
    echo -e "  -p  Number of subprocesses for multiprocessing tools.\n" \
            "     Defaults to the available cores."
    echo -e "  -i  Specify to allow partial alignments."
    echo -e "  -z  Specify to disable cleaning and compression of temporary files."
    echo -e "  -h  This message."
}

# Argument validation variables
READS_SET=false
A10X_SET=false
REFDIR_SET=false
SPECIES_SET=false
LOCI_SET=false
DIST_SET=false
MODEL_SET=false
THRESHOLD_METHOD_SET=false
THRESHOLD_MODEL_SET=false
CUTOFF_SET=false
IGDATA_SET=false
OUTNAME_SET=false
OUTDIR_SET=false
FORMAT_SET=false
NPROC_SET=false

# Argument defaults
PARTIAL=""
ZIP_FILES=true
DELETE_FILES=true
SPC=0.995

# Get commandline arguments
while getopts "s:a:r:g:t:x:m:e:d:u:b:n:o:f:p:izh" OPT; do
    case "$OPT" in
    s)  READS=$OPTARG
        READS_SET=true
        ;;
    a)  A10X=$OPTARG
        A10X_SET=true
        ;;
    r)  REFDIR=$OPTARG
        REFDIR_SET=true
        ;;
    g)  SPECIES=$OPTARG
        SPECIES_SET=true
        ;;
    t)  LOCI=$OPTARG
        LOCI_SET=true
        ;;
    x)  DIST=$OPTARG
        DIST_SET=true
        ;;
    m)  MODEL=$OPTARG
        MODEL_SET=true
        ;;
    e)  THRESHOLD_METHOD=$OPTARG
        THRESHOLD_METHOD_SET=true
        ;;
    d)  THRESHOLD_MODEL=$OPTARG
        THRESHOLD_MODEL_SET=true
        ;;
    u)  CUTOFF=$OPTARG
        CUTOFF_SET=true
        ;;
    b)  IGDATA=$OPTARG
        IGDATA_SET=true
        ;;
    n)  OUTNAME=$OPTARG
        OUTNAME_SET=true
        ;;
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    f)  FORMAT=$OPTARG
        FORMAT_SET=true
        ;;
    p)  NPROC=$OPTARG
        NPROC_SET=true
        ;;
    i)  PARTIAL="--partial"
        ;;
    z)  ZIP_FILES=false
        DELETE_FILES=false
        ;;
    h)  print_usage
        exit
        ;;
    \?) echo -e "Invalid option: -${OPTARG}" >&2
        exit 1
        ;;
    :)  echo -e "Option -${OPTARG} requires an argument" >&2
        exit 1
        ;;
    esac
done

# Exit if required arguments are not provided
if ! ${READS_SET}; then
    echo -e "You must specify the input sequences using the -s option." >&2
    exit 1
fi

if ! ${A10X_SET}; then
    echo -e "You must specify the Cell Ranger annotation file using the -x option." >&2
    exit 1
fi

# Check that files exist and determined absolute paths
if [ -e ${READS} ]; then
    READS=$(realpath ${READS})
else
    echo -e "File '${READS}' not found." >&2
    exit 1
fi

if [ -e ${A10X} ]; then
    A10X=$(realpath ${A10X})
else
    echo -e "File '${A10X}' not found." >&2
    exit 1
fi

# Set and check species
if ! ${SPECIES_SET}; then
    SPECIES="human"
elif [ ${SPECIES} != "human" ] && \
     [ ${SPECIES} != "mouse" ] && \
     [ ${SPECIES} != "rabbit" ] && \
     [ ${SPECIES} != "rat" ] && \
     [ ${SPECIES} != "rhesus_monkey" ]; then
    echo "Species (-g) must be one of 'human', 'mouse', 'rabbit', 'rat', or 'rhesus_monkey'." >&2
    exit 1
fi

# Set regions
REGIONS="default"

# Set and check receptor type
if ! ${LOCI_SET}; then
    LOCI="ig"
elif [ ${LOCI} != "ig" ] && [ ${LOCI} != "tr" ]; then
    echo "Receptor type (-t) must be one of 'ig' or 'tr'." >&2
    exit 1
fi

# Set reference sequence
if ! ${REFDIR_SET}; then
    REFDIR="/usr/local/share/germlines/imgt/${SPECIES}/vdj"
else
    REFDIR=$(realpath ${REFDIR})
fi

# Set distance model
if ! ${MODEL_SET}; then
    MODEL="nt"
fi

# Set threshold method
if ! ${THRESHOLD_METHOD_SET}; then
    THRESHOLD_METHOD="density"
fi

# Set threshold model
if ! ${THRESHOLD_MODEL_SET}; then
    THRESHOLD_MODEL="gamma-gamma"
fi


# Set cutoff (method to use for threshold selection)
if ! ${CUTOFF_SET}; then
    CUTOFF="user"
fi

# Set blast database
if ! ${IGDATA_SET}; then
    IGDATA="/usr/local/share/igblast"
else
    IGDATA=$(realpath ${IGDATA})
fi

# Set output name
if ! ${OUTNAME_SET}; then
    OUTNAME=$(basename ${READS} | sed 's/\.[^.]*$//; s/_L[0-9]*_R[0-9]_[0-9]*//')
fi

# Set output directory
if ! ${OUTDIR_SET}; then
    OUTDIR=${OUTNAME}
fi

# Check output directory permissions
if [ -e ${OUTDIR} ]; then
    if ! [ -w ${OUTDIR} ]; then
        echo -e "Output directory '${OUTDIR}' is not writable." >&2
        exit 1
    fi
else
    PARENTDIR=$(dirname $(realpath ${OUTDIR}))
    if ! [ -w ${PARENTDIR} ]; then
        echo -e "Parent directory '${PARENTDIR}' of new output directory '${OUTDIR}' is not writable." >&2
        exit 1
    fi
fi


# Set format options
if ! ${FORMAT_SET}; then
    FORMAT="airr"
fi

if [[ "${FORMAT}" == "airr" ]]; then
    EXT="tsv"
    LOCUS_FIELD="locus"
    PROD_FIELD="productive"
else
	EXT="tab"
	LOCUS_FIELD="LOCUS"
	PROD_FIELD="FUNCTIONAL"
fi

# Set number of processes
if ! ${NPROC_SET}; then
    NPROC=$(nproc)
fi

# Define pipeline steps
SPLIT=true
GERMLINES=true
if ! ${DIST_SET}; then
    CLONE=false
else
    CLONE=true
fi

# DefineClones run parameters
DC_MODE="gene"
DC_ACT="set"

# Create germlines parameters
CG_GERM="full dmask"

# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="logs"
PIPELINE_LOG="${LOGDIR}/pipeline-10x.log"
ERROR_LOG="${LOGDIR}/pipeline-10x.err"
mkdir -p ${LOGDIR}
echo '' > $PIPELINE_LOG
echo '' > $ERROR_LOG

# Check for errors
check_error() {
    if [ -s $ERROR_LOG ]; then
        echo -e "ERROR:"
        cat $ERROR_LOG | sed 's/^/    /'
        exit 1
    fi
}

# Set extension
IGBLAST_VERSION=$(igblastn -version  | grep 'Package' |sed s/'Package: '//)
CHANGEO_VERSION=$(python3 -c "import changeo; print('%s-%s' % (changeo.__version__, changeo.__date__))")
SCOPER_VERSION=$(Rscript -e "v <- packageVersion('scoper');cat(as.character(v))")

# Start
echo -e "IDENTIFIER: ${OUTNAME}"
echo -e "DIRECTORY: ${OUTDIR}"
echo -e "CHANGEO VERSION: ${CHANGEO_VERSION}"
echo -e "IGBLAST VERSION: ${IGBLAST_VERSION}"
echo -e "SCOPER_VERSION: ${SCOPER_VERSION}"
echo -e "\nSTART"
STEP=0

# Convert to FASTA if needed
BASE_NAME=$(basename ${READS})
EXT_NAME=${BASE_NAME##*.}
if [ "${EXT_NAME,,}" == "fastq" ] || [ "${EXT_NAME,,}" == "fq" ]; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Convert to FASTA"
    IG_FILE=$(fastq2fasta.py ${READS})
else
    IG_FILE=${READS}
fi

# Run IgBLAST
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "AssignGenes igblast"
AssignGenes.py igblast -s ${IG_FILE} --organism ${SPECIES} --loci ${LOCI} \
    -b ${IGDATA} --format blast --nproc ${NPROC} \
    --outname "${OUTNAME}" --outdir . \
     >> $PIPELINE_LOG 2> $ERROR_LOG
FMT7_FILE="${OUTNAME}_igblast.fmt7"
check_error

# Parse IgBLAST output
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "MakeDb igblast"
MakeDb.py igblast -i ${FMT7_FILE} -s ${IG_FILE} --10x ${A10X} -r ${REFDIR} \
    --extended --failed ${PARTIAL} --outname "${OUTNAME}" --format ${FORMAT} \
    --regions ${REGIONS} \
    >> $PIPELINE_LOG 2> $ERROR_LOG
DB_PASS="${OUTNAME}_db-pass.${EXT}"
DB_FAIL="${OUTNAME}_db-fail.${EXT}"
check_error

# Split by chain and productivity
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "ParseDb select"
ParseDb.py select -d ${DB_PASS} -f ${LOCUS_FIELD} -u IGH TRB TRD \
    -o "${OUTNAME}_heavy.${EXT}" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
ParseDb.py select -d ${DB_PASS} -f ${LOCUS_FIELD} -u IGK IGL TRA TRG \
    -o "${OUTNAME}_light.${EXT}" \
    >> $PIPELINE_LOG 2> $ERROR_LOG
HEAVY_ALL="${OUTNAME}_heavy.${EXT}"
LIGHT_ALL="${OUTNAME}_light.${EXT}"

printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "ParseDb split"
ParseDb.py split -d "${OUTNAME}_heavy.${EXT}" -f ${PROD_FIELD} \
    >> $PIPELINE_LOG 2> $ERROR_LOG
ParseDb.py split -d "${OUTNAME}_light.${EXT}" -f ${PROD_FIELD} \
    >> $PIPELINE_LOG 2> $ERROR_LOG
HEAVY_PROD="${OUTNAME}_heavy_${PROD_FIELD}-T.${EXT}"
LIGHT_PROD="${OUTNAME}_light_${PROD_FIELD}-T.${EXT}"
HEAVY_NON="${OUTNAME}_heavy_${PROD_FIELD}-F.${EXT}"
LIGHT_NON="${OUTNAME}_light_${PROD_FIELD}-F.${EXT}"

# Assign clones
if $CLONE; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Single cell filter"
    singlecell-filter -d ${HEAVY_PROD},${LIGHT_PROD} -o . -f ${FORMAT} \
    >> $PIPELINE_LOG 2> $ERROR_LOG
    check_error

    HEAVY_PROD="${OUTNAME}_heavy_${PROD_FIELD}-T_sc-pass.${EXT}"
    LIGHT_PROD="${OUTNAME}_light_${PROD_FIELD}-T_sc-pass.${EXT}"
    if [ "$DIST" == "auto" ]; then
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Detect cloning threshold"
        shazam-threshold -d ${HEAVY_PROD},${LIGHT_PROD}  -m ${THRESHOLD_METHOD} -n "${OUTNAME}" \
        --model ${THRESHOLD_MODEL} --cutoff ${CUTOFF} --spc ${SPC} -o . \
        -f ${FORMAT} -p ${NPROC} \
        > /dev/null 2> $ERROR_LOG
        check_error
        DIST=$(tail -n1 "${OUTNAME}_threshold-values.tab" | cut -f2)
    else
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Calculating distances"
        shazam-threshold -d ${HEAVY_PROD} -m none -n "${OUTNAME}" -o . \
        -f ${FORMAT} -p ${NPROC} \
        > /dev/null 2> $ERROR_LOG
        check_error
    fi

    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Define clones (scoper)"
    scoper-clone -d ${HEAVY_PROD},${LIGHT_PROD} -o . -f ${FORMAT} \
        --method ${MODEL} --threshold ${DIST} --nproc ${NPROC} \
        --log "${LOGDIR}/clone.log" \
        --name "${OUTNAME}_heavy","${OUTNAME}_light" \
        >> $PIPELINE_LOG 2> $ERROR_LOG

    CLONE_FILE="${OUTNAME}_heavy_clone-pass.${EXT}"
    check_error

    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "CreateGermlines"
    CreateGermlines.py -d ${CLONE_FILE} --cloned -r ${REFDIR} -g ${CG_GERM} \
        --outname "${OUTNAME}_heavy" --log "${LOGDIR}/germline.log" --format ${FORMAT} \
        >> $PIPELINE_LOG 2> $ERROR_LOG
    HEAVY_PROD="${OUTNAME}_heavy_germ-pass.${EXT}"
	check_error
fi

# Zip or delete intermediate files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Compressing files"
TEMP_FILES=$(ls *.tsv *.tab 2> /dev/null | grep -v "${HEAVY_PROD}\|${LIGHT_PROD}\|${HEAVY_NON}\|${LIGHT_NON}\|${DB_PASS}")
if [[ ! -z $TEMP_FILES ]]; then
    if $ZIP_FILES; then
        tar -zcf temp_files.tar.gz $TEMP_FILES
    fi
    if $DELETE_FILES; then
        rm $TEMP_FILES
    fi
fi

# End
printf "DONE\n\n"
cd ../
