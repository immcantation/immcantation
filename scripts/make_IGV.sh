#!/usr/bin/env bash
# 
# make gapped and ungapped IGV germline fasta files
# from the IGHV and IGKV fasta files in the 
# VDJ folder
# 
# Author:  Susanna Marquez
# Date:    2016.05.02
#
# Arguments:
#   -v = Path to dir with IMGT IGHV and IGKV fasta files
#   -o = Output directory.
#   -h = Display help.

# Default argument values
OUTDIR="."

# Print usage
usage () {
    echo "Usage: `basename $0` [OPTIONS]"
    echo "  -v  Path to dir with IMGT IGHV and IGKV fasta files."
    echo "  -o  Output directory."
    echo "  -h  This message."
}

# Get commandline arguments
while getopts "v:o:h" OPT; do
    case "$OPT" in
    v)  VDIR=$(readlink -f $OPTARG)
        VDIR_SET=true
        ;;
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    h)  usage
        exit
        ;;
    \?) echo "Invalid option $OPTARG" >&2
        exit 1
        ;;
    :)  echo "Option $OPTARG requires an argument" >&2
        exit 1
        ;;
    esac
done

# Exit if dir doesn't exist
if [ ! -d "$VDIR" ]; then
   echo "$VDIR doesn't exist" >&2
   exit 1
fi

IGHV=$(find $VDIR -type f -name  *IGHV.fasta)
NUMFILES=${#IGHV[@]}

if [ "$NUMFILES" -eq "0" ]; then
	# Exit if no expected input files found in dir
	echo "No *IGHV.fasta files found\n"
	exit 1
elif [ "$NUMFILES" -gt "1" ]; then
	# Exit if more than one expected input files found in dir
	echo "More than one *IGHV.fasta files found\n"
	exit 1
fi

IGKV=$(find $VDIR -type f -name  *IGKV.fasta)
NUMFILES=${#IGKV[@]}

if [ "$NUMFILES" -eq "0" ]; then
	# Exit if no expected input files found in dir
	echo "No *IGKV.fasta files found\n"
	exit 1
elif [ "$NUMFILES" -gt "1" ]; then
	# Exit if more than one expected input files found in dir
	echo "More than one *IGKV.fasta files found\n"
	exit 1
fi

# Make output directory if it does not exist
if $OUTDIR_SET && [ ! -d "${OUTDIR}" ]; then
     mkdir -p $OUTDIR
fi

OUTFILE=$(basename $IGHV)

GAPPED=${OUTFILE//_IGHV.fasta/IGV_gapped.fasta}
GAPPED=$OUTDIR"/"$GAPPED

UNGAPPED=${OUTFILE//_IGHV.fasta/IGV_ungappped.fasta}
UNGAPPED=$OUTDIR"/"$UNGAPPED

## Concatenate
cat $IGHV $IGKV > $GAPPED

## Ungap
command -v seqmagick >/dev/null 2>&1 || { echo >&2 "Can't find seqmagick."; exit 1; }
seqmagick convert --ungap $GAPPED $UNGAPPED
