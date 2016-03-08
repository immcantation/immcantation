#!/usr/bin/env bash
# Simple IgBLAST wrapper
#
# Author:  Jason Anthony Vander Heiden
# Date:    2016.03.06
#
# Arguments:
#   -s = FASTA sequence file.
#   -o = Output directory.
#   -n = Number of IgBLAST threads.
#   -h = Display help.

# Default argument values
OUTDIR="."
NPROC=1

# Define IgBLAST directories and base command
export IGDATA="${HOME}/apps/igblast-1.4.0"
IGBLAST_DB="${IGDATA}/database"
IGBLAST_CMD="${IGDATA}/igblastn \
    -germline_db_V ${IGBLAST_DB}/imgt_human_IG_V \
    -germline_db_D ${IGBLAST_DB}/imgt_human_IG_D \
    -germline_db_J ${IGBLAST_DB}/imgt_human_IG_J \
    -auxiliary_data ${IGDATA}/optional_file/human_gl.aux \
    -domain_system imgt -ig_seqtype Ig -organism human \
    -outfmt '7 std qseq sseq btop'"

# Print usage
usage () {
    echo "Usage: `basename $0` [OPTIONS]"
    echo "  -s  FASTA sequence file."
    echo "  -o  Output directory."
    echo "  -n  Number of IgBLAST threads."
    echo "  -h  This message."
}

# Get commandline arguments
while getopts "s:o:n:e:h" OPT; do
    case "$OPT" in
    s)  READFILE=$(readlink -f $OPTARG)
        READFILE_SET=true
        ;;
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    n)  NPROC=$OPTARG
        ;;
    h)  usage
        exit
        ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)  echo "Option -$OPTARG requires an argument" >&2
        exit 1
        ;;
    esac
done

# Exit if no input file provided
if ! $READFILE_SET; then
    echo "You must specify a FASTA file using the -s option" >&2
    exit 1
fi

# Make output directory if it does not exist
if $OUTDIR_SET %% [ ! -d "${OUTDIR}" ]; then
    mkdir -p $OUTDIR
fi

# Set run commmand
OUTFILE=$(basename ${READFILE})
OUTFILE="${OUTDIR}/${OUTFILE%.fasta}.fmt7"
IGBLAST_VER=$(${IGBLAST_CMD} -version | grep 'Package' |sed s/'Package: '//)
IGBLAST_RUN="${IGBLAST_CMD} -query ${READFILE} -out ${OUTFILE} -num_threads ${NPROC}"

# Align V(D)J segments using IgBLAST
echo -e "   START> igblastn"
echo -e " VERSION> ${IGBLAST_VER}"
echo -e "  GERMDB> ${IGBLAST_DB}"
echo -e "    FILE> $(basename ${READFILE})\n"
echo -e "PROGRESS> [Running]"
eval $IGBLAST_RUN
echo -e "PROGRESS> [Done   ]\n"
echo -e "  OUTPUT> $(basename ${OUTFILE})"
echo -e "     END> igblastn\n"