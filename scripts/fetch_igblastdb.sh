#!/usr/bin/env bash
# Download IgBLAST database files
#
# Author:  Jason Vander Heiden
# Date:    2020.08.11
#
# Arguments:
#   -o = Output directory for downloaded files. Defaults to current directory.
#   -x = Exclude download of internal_data and optional_file bundles.
#   -h = Display help.

# Default argument values
OUTDIR="."
DOWNLOAD_ALL=true

# Print usage
usage () {
    echo "Usage: `basename $0` [OPTIONS]"
    echo "  -o  Output directory for downloaded files. Defaults to current directory."
    echo "  -x  Exclude download of internal_data and optional_file bundles."
    echo "  -h  This message."
}

# Get commandline arguments
while getopts "o:xh" OPT; do
    case "$OPT" in
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    x)  DOWNLOAD_ALL=false
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

# Fetch database
wget -q -r -nH --cut-dirs=5 --no-parent \
    ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/database \
    -P ${OUTDIR}/database

if $DOWNLOAD_ALL; then
    # Fetch internal_data
    wget -q -r -nH --cut-dirs=5 --no-parent \
        ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data \
        -P ${OUTDIR}/internal_data

    # Fetch optional_file
    wget -q -r -nH --cut-dirs=5 --no-parent \
        ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/optional_file  \
        -P ${OUTDIR}/optional_file
fi