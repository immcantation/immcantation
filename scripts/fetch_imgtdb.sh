#!/usr/bin/env bash
# Download germlines from the IMGT website
#
# Author:  Mohamed Uduman, Jason Anthony Vander Heiden
# Date:    2017.04.06
#
# Arguments:
#   -o = Output directory for downloaded files. Defaults to current directory.
#   -h = Display help.

# Default argument values
OUTDIR="."

# Print usage
usage () {
    echo "Usage: `basename $0` [OPTIONS]"
    echo "  -o  Output directory for downloaded files. Defaults to current directory."
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

# IMGT Gapped Germlines
REPERTOIRE="imgt"

# Associative array where keys are species folder names and values are query strings
declare -A SPECIES_QUERY
SPECIES_QUERY[human]="Homo+sapiens"
SPECIES_QUERY[mouse]="Mus"

# Associative a
declare -A SPECIES_REPLACE
SPECIES_REPLACE[human]='s/Homo sapiens/Homo_sapiens/g'
SPECIES_REPLACE[mouse]='s/Mus musculus/Mus_musculus/g'

# For each species
for KEY in ${!SPECIES_QUERY[@]}
do
	echo "Downloading ${KEY} repertoires into ${OUTDIR}..."

	# Download VDJ
	echo "|- VDJ regions"
    FILE_PATH="${OUTDIR}/${KEY}/vdj"
    mkdir -p $FILE_PATH

    # VDJ Ig
    echo "|---- Ig"
    for CHAIN in IGHV IGHD IGHJ IGKV IGKJ IGLV IGLJ
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.14+${CHAIN}&species=${SPECIES_QUERY[${KEY}]}"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i "${SPECIES_REPLACE[${KEY}]}" $FILE_NAME
        rm $TMP_FILE
    done

    # VDJ TCR
    echo "|---- TCR"
    for CHAIN in TRAV TRAJ TRBV TRBD TRBJ TRDV TRDD TRDJ TRGV TRGJ
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.14+${CHAIN}&species=${SPECIES_QUERY[${KEY}]}"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i "${SPECIES_REPLACE[${KEY}]}" $FILE_NAME
        rm $TMP_FILE
    done


	# Download leaders
    echo "|- Spliced leader regions"
    FILE_PATH="${OUTDIR}/${KEY}/leader"
    mkdir -p $FILE_PATH

    # Leader Ig
    echo "|---- Ig"
    for CHAIN in IGH IGK IGL
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=8.1+${CHAIN}V&species=${SPECIES_QUERY[${KEY}]}&IMGTlabel=L-PART1+L-PART2"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}L.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i "${SPECIES_REPLACE[${KEY}]}" $FILE_NAME
        rm $TMP_FILE
    done

    # Leader TCR
    echo "|---- TCR"
    for CHAIN in TRA TRB TRG TRD
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=8.1+${CHAIN}V&species=${SPECIES_QUERY[${KEY}]}&IMGTlabel=L-PART1+L-PART2"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}L.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i "${SPECIES_REPLACE[${KEY}]}" $FILE_NAME
        rm $TMP_FILE
    done


	# Download constant regions
    echo "|- Spliced constant regions"
    FILE_PATH="${OUTDIR}/${KEY}/constant/"
    mkdir -p $FILE_PATH

    # Constant Ig
    echo "|---- Ig"
    for CHAIN in IGHC IGKC IGLC
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=14.1+${CHAIN}&species=${SPECIES_QUERY[${KEY}]}"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i "${SPECIES_REPLACE[${KEY}]}" $FILE_NAME
        rm $TMP_FILE
    done

    # Constant for TCR
    echo "|---- TCR"
    for CHAIN in TRAC TRBC TRGC TRDC
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=14.1+${CHAIN}&species=${SPECIES_QUERY[${KEY}]}"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i "${SPECIES_REPLACE[${KEY}]}" $FILE_NAME
        rm $TMP_FILE
    done

    echo ""

done
