#!/usr/bin/env bash
# Super script to preprocess fastq files for pRESTO.
# Will blast reads against Phi-X174 genome. Reads with 
# a hit will be filtered out
#
# Author:  Susanna Marquez
# Date:    2018.03.19
#
# Arguments:
#   -s  FASTQ sequence file.
#   -r  Directory containing phiX174 reference db.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the input filename.
#   -o  Output directory. Will be created if it does not exist.
#       Defaults to a directory matching the sample identifier in the current working directory.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -h  Display help

# Print usage
print_usage() {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -s   FASTQ sequence file."
    echo -e "  -r   Directory containing phiX174 reference db.\n" \
            "      Defaults to /usr/local/share/phix."
    echo -e "  -n   Sample identifier which will be used as the output file prefix.\n" \
            "      Defaults to a truncated version of the input filename."
    echo -e "  -o  Output directory. Will be created if it does not exist.\n" \
            "     Defaults to a directory matching the sample identifier in the current working directory."
    echo -e "  -p   Number of subprocesses for multiprocessing tools.\n" \
            "      Defaults to the available cores."
    echo -e "  -h   This message."
}


# Argument validation variables
READS_SET=false
PHIXDIR_SET=false
OUTDIR_SET=false
OUTNAME_SET=false
NPROC_SET=false

# Define BLAST command
BLAST="blastn"

# Get commandline arguments
while getopts "s:r:n:o:p:h" OPT; do
    case "$OPT" in
    s)  READS=$OPTARG
        READS_SET=true
        ;;
    r)  PHIXDIR=$OPTARG
        PHIXDIR_SET=true
        ;;
    n)  OUTNAME=$OPTARG
        OUTNAME_SET=true
        ;;
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    p)  NPROC=$OPTARG
        NPROC_SET=true
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

# Check that files exist and determined absolute paths
if [ -e ${READS} ]; then
    READS=$(realpath ${READS})
else
    echo -e "File '${READS}' not found." >&2
    exit 1
fi

# Exit if required arguments are not provided
if ! ${PHIXDIR_SET}; then
    PHIXDIR="/usr/local/share/phix"
fi

# Check that dir exists and determined absolute paths
if [ -e ${PHIXDIR} ]; then
    PHIXDIR=$(realpath ${PHIXDIR})
    PHIXDB=$(ls ${PHIXDIR}/*fna)
else
    echo -e "Directory '${PHIXDIR}' not found." >&2
    exit 1
fi

# Set output name
if ! ${OUTNAME_SET}; then
     #OUTNAME=$(basename ${READS} | sed 's/.fastq/_nophix/')
     OUTNAME=$(basename ${READS} | sed 's/.fastq//')
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

# Set number of processes
if ! ${NPROC_SET}; then
    NPROC=$(nproc)
fi

# Make output directory
mkdir -p ${OUTDIR}; cd ${OUTDIR}

# Define log files
LOGDIR="logs"
PIPELINE_LOG="${LOGDIR}/pipeline-phix.log"
ERROR_LOG="${LOGDIR}/pipeline-phix.err"
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

# Start
BLASTN_VERSION=$(blastn -version  | grep 'Package' |sed s/'Package: '//)
PHIX_VERSION=$(grep date ${PHIXDIR}/PhiX174.yaml | sed s/'date: *'//)

echo -e "OUTNAME: ${OUTNAME}"
echo -e "OUTDIR: ${OUTDIR}"
echo -e "PHIXDB:  ${PHIXDB}"
echo -e "BLASTN VERSION: ${BLASTN_VERSION}"
echo -e "PHIX VERSION (DOWNLOAD DATE): ${PHIX_VERSION}"
echo -e "LOGDIR: ${LOGDIR}"
echo -e "\nSTART"
STEP=0

# Remove all-N sequence becase blastn crashes with all N sequences)
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Removing all N sequences"
echo -e "       START> awk" >> $PIPELINE_LOG
NO_N_READS="${OUTNAME}_clean-missing.fastq"
awk '{y= i++ % 4 ; L[y]=$0; if(y==3 && L[1] ~ /[^N]/) {printf("%s\n%s\n%s\n%s\n",L[0],L[1],L[2],L[3]);}}' ${READS} \
    > ${NO_N_READS} 2> $ERROR_LOG

# Set input file for alignment
INPUT_SIZE=$((`wc -l < ${READS}`/4))
OUTPUT_SIZE=$((`wc -l < ${NO_N_READS}`/4))
REMOVED_SEQS=$((${INPUT_SIZE}-${OUTPUT_SIZE}))

if [ ${REMOVED_SEQS} -gt 0 ]; then
   CONVERT_FILE=${NO_N_READS}
else
   CONVERT_FILE=${READS}
   rm $NO_N_READS
fi
   
echo -e "  INPUT_SIZE> ${INPUT_SIZE}" >> $PIPELINE_LOG
echo -e " OUTPUT_SIZE> ${OUTPUT_SIZE}" >> $PIPELINE_LOG
echo -e "REMOVED_SEQS> ${REMOVED_SEQS}" >> $PIPELINE_LOG
echo -e "  READS_FILE> ${READS}\n" >> $PIPELINE_LOG

# Convert headers to presto format
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ConvertHeaders"
ConvertHeaders.py illumina -s ${CONVERT_FILE} --outdir ${OUTDIR} --outname ${OUTNAME} --fasta \
    >> $PIPELINE_LOG 2> $ERROR_LOG
FASTA_FILE="${OUTNAME}_convert-pass.fasta"
check_error

# Run blastn
BLAST_CMD="${BLAST} \
     -query ${FASTA_FILE} \
     -db ${PHIXDB} \
     -outfmt '6 std qseq sseq btop' \
     -out ${OUTNAME}_phix.fmt6 \
     -num_threads ${NPROC}"

printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "BLASTN"
echo -e "   START> blastn" >> $PIPELINE_LOG
echo -e "    FILE> $(basename ${FASTA_FILE}) \n" >> $PIPELINE_LOG
echo -e "PROGRESS> [Running]" >> $PIPELINE_LOG
eval ${BLAST_CMD} >> $PIPELINE_LOG 2> $ERROR_LOG
echo -e "PROGRESS> [Done   ]\n" >> $PIPELINE_LOG
echo -e "  OUTPUT> ${OUTNAME}_phix.fmt6" >> $PIPELINE_LOG
echo -e "     END> blastn\n" >> $PIPELINE_LOG
check_error

# Add header, need ID column name for Splitseq
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Add header"
sed -i '1iID' "${OUTNAME}_phix.fmt6"
ID_FILE="${OUTNAME}_phixhits.txt"
sed -r '2,$ s/(^[^\|]*).*/\1/' "${OUTNAME}_phix.fmt6" > ${ID_FILE}

# Filter input fasta/q to names not in the .fmt6 file
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq select"
SplitSeq.py select -s ${READS} -f ID -t ${ID_FILE} --not --outdir ${OUTDIR} --outname ${OUTNAME} \
    >> $PIPELINE_LOG 2> $ERROR_LOG
check_error

# Remove temporary files
rm $FASTA_FILE
if [ ${REMOVED_SEQS} -gt 0 ]; then
   rm $NO_N_READS
fi

# End
printf "DONE\n\n"
