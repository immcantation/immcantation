#!/bin/bash

IGBLAST_PATH=~/apps/igblast-1.4.0
GERMLINES_PATH=~/workspace/igpipeline/germlines/IMGT/Human/VDJ

cd $IGBLAST_PATH

# Create fasta files of each chain/segment combo
for CHAIN in IG TR
do
	for SEGMENT in V D J
	do
		cat ${GERMLINES_PATH}/IMGT_Human_${CHAIN}?${SEGMENT}.fasta > ${IGBLAST_PATH}/imgt_human_${CHAIN}_${SEGMENT}.fasta
	done
done

# Parse each created fasta file to create igblast database
for F in $(ls *.fasta)
do
	./edit_imgt_file.pl $F > database/${F%%.*}
	./makeblastdb -parse_seqids -dbtype nucl -in database/${F%%.*}
	rm $F
done
