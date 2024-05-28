#!/usr/bin/env bash
# Download IgBLAST database files
#
# Author:  Jason Vander Heiden
# Date:    2020.08.11
#
# Arguments:
#   -o = Output directory for downloaded files. Defaults to current directory.
#   -x = Include download of internal_data and optional_file bundles.
#   -h = Display help.

# Default argument values
OUTDIR="."
DOWNLOAD_ALL=false

# Print usage
usage () {
    echo "Usage: `basename $0` [OPTIONS]"
    echo "  -o  Output directory for downloaded files. Defaults to current directory."
    echo "  -x  Include download of legacy internal_data and optional_file bundles."
    echo "  -h  This message."
}

# Get commandline arguments
while getopts "o:xh" OPT; do
    case "$OPT" in
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    x)  DOWNLOAD_ALL=true
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

# Make output directory if it does not exist
if $OUTDIR_SET && [ ! -d "${OUTDIR}" ]; then
    mkdir -p $OUTDIR
fi

WGET="wget"
if ! command -v wget &> /dev/null; then
    if ! command -v wget2 &> /dev/null; then
        echo "wget or wget2 not found."
        exit 1
    else
        WGET=wget2
    fi
fi
WGET_MAJOR=$(${WGET} --version | grep "Wget[0-9]* [0-9.]*" | sed 's/^[^0-9]*\([0-9]*\).*$/\1/')
if [ "$WGET_MAJOR" -lt 2 ]; then
    ROBOTS=("--execute robots=off")
else
    ROBOTS="--no-robots"
fi

# Fetch database
${WGET} -q -r -nH --cut-dirs=5 --no-parent ${ROBOTS} \
        https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/database/ \
        -P ${OUTDIR}/database

# Check one file to verify download worked 
if [ ! -s "${OUTDIR}/database/mouse_gl_VDJ.tar" ]
then
    echo "${OUTDIR}/database/mouse_gl_VDJ.tar. Is ncbi.nlm.nih.gov server running?"
    exit 1
fi

# Extract
tar -C ${OUTDIR}/database -xf ${OUTDIR}/database/mouse_gl_VDJ.tar
tar -C ${OUTDIR}/database -xf ${OUTDIR}/database/rhesus_monkey_VJ.tar

if $DOWNLOAD_ALL; then
    # Fetch internal_data
    ${WGET} -q -r -nH --cut-dirs=5 --no-parent ${ROBOTS} \
        https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/old_internal_data/ \
        -P ${OUTDIR}/internal_data
    # Check one file to verify download worked 
    if [ ! -s "${OUTDIR}/internal_data/human/human_V.nhr" ]
    then
        echo "${OUTDIR}/internal_data/human/human_V.nhr not found. Is ncbi.nlm.nih.gov server running?"
        exit 1
    fi

    # Fetch optional_file
    ${WGET} -q -r -nH --cut-dirs=5 --no-parent ${ROBOTS} \
        https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/old_optional_file/ \
        -P ${OUTDIR}/optional_file
    # Check one file to verify download worked 
    if [ ! -s "${OUTDIR}/optional_file/human_gl.aux" ]
    then
        echo "${OUTDIR}/optional_file/human_gl.aux not found. Is ncbi.nlm.nih.gov server running?"
        exit 1
    fi        
fi