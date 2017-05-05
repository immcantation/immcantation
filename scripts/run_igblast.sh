#!/usr/bin/env bash
# Simple IgBLAST wrapper
#
# Author:  Jason Anthony Vander Heiden
# Date:    2016.05.05
#
# Arguments:
#   -s = FASTA sequence file.
#   -o = Output directory.
#   -d = IGDATA directory, which contains the IgBLAST database, optional_file
#        and auxillary_data directories. Defaults to /usr/local/share/igblast.
#   -n = Number of IgBLAST threads. Defaults to 1.
#   -h = Display help.

# Default argument values
IGDATA="/usr/local/share/igblast"
OUTDIR="."
NPROC=1

# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -s  FASTA sequence file."
    echo -e "  -o  Output directory."
    echo -e "  -d  IGDATA directory, which contains the IgBLAST database,\n" \
            "     optional_file and auxillary_data directories.\n" \
            "     Defaults to /usr/local/share/igblast."
    echo -e "  -n  Number of IgBLAST threads. Defaults to 1."
    echo -e "  -h  This message."
}

# Validation variables
READFILE_SET=false
OUTDIR_SET=false

# Get commandline arguments
while getopts "s:o:d:n:h" OPT; do
    case "$OPT" in
    s)  READFILE=$(readlink -f $OPTARG)
        READFILE_SET=true
        ;;
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    d)  IGDATA=$OPTARG
        IGDATA_SET=true
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
if $OUTDIR_SET && [ ! -d "${OUTDIR}" ]; then
    mkdir -p $OUTDIR
fi

# Define IgBLAST directories and base command
export IGDATA
IGBLAST_DB="${IGDATA}/database"
IGBLAST_CMD="igblastn \
    -germline_db_V ${IGBLAST_DB}/imgt_human_ig_v \
    -germline_db_D ${IGBLAST_DB}/imgt_human_ig_d \
    -germline_db_J ${IGBLAST_DB}/imgt_human_ig_j \
    -auxiliary_data ${IGDATA}/optional_file/human_gl.aux \
    -domain_system imgt -ig_seqtype Ig -organism human \
    -outfmt '7 std qseq sseq btop'"

# Set run commmand
OUTFILE=$(basename ${READFILE})
OUTFILE="${OUTDIR}/${OUTFILE%.fasta}.fmt7"
IGBLAST_VER=$(${IGBLAST_CMD} -version | grep 'Package' |sed s/'Package: '//)
IGBLAST_RUN="${IGBLAST_CMD} -query ${READFILE} -out ${OUTFILE} -num_threads ${NPROC}"

# Align V(D)J segments using IgBLAST
echo -e "   START> igblastn"
echo -e " VERSION> ${IGBLAST_VER}"
echo -e "  IGDATA> ${IGDATA}"
echo -e "  GERMDB> imgt_human"
echo -e "    FILE> $(basename ${READFILE})\n"
echo -e "PROGRESS> [Running]"
eval $IGBLAST_RUN
echo -e "PROGRESS> [Done   ]\n"
echo -e "  OUTPUT> $(basename ${OUTFILE})"
echo -e "     END> igblastn\n"