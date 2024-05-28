#!/usr/bin/env bash
# Download PhiX174 nucleotide sequence from NCBI and build blastn indexed db
#
# Author:  Susanna Marquez
# Date:    2017.09.21
#
# Arguments:
#   -o = Output directory. Defaults to current directory.
#   -h = Display help.

# Default argument values
OUTDIR="."

# Print usage
usage () {
    echo "Usage: `basename $0` [OPTIONS]"
    echo "  -o  Output directory. Defaults to current directory."
    echo "  -h  This message."
}

# Get commandline arguments
while getopts "o:h" OPT; do
    case "$OPT" in
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
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

# Info
DATE=$(date +"%Y.%m.%d")

WGET="wget"
if ! command -v wget &> /dev/null; then
    if ! command -v wget2 &> /dev/null; then
        echo "wget or wget2 not found."
        exit 1
    else
        WGET=wget2
    fi
fi

# Download and unpack
${WGET} "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz" -P $OUTDIR
gunzip "${OUTDIR}/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
BLAST_DB="${OUTDIR}/GCF_000819615.1_ViralProj14015_genomic.fna"

# Build blast database
makeblastdb -in ${BLAST_DB} -parse_seqids -dbtype nucl

# Write download info
INFO_FILE="${OUTDIR}/PhiX174.yaml"
echo -e "source:  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz" > $INFO_FILE
echo -e "date:    ${DATE}" >> $INFO_FILE
