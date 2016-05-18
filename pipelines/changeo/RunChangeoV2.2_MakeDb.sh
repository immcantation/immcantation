#!/bin/bash
# Script to run MakeDb and ParseDb on multiple inputs
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2015.07.08
#
# Required Arguments:
#   $1 = A two column tab delimited file mapping IMGT zip files (column 1) 
#        to submitted FASTA files (column 2)
#   $2 = Output directory

mkdir -p $2/db
mkdir -p $2/functional

while read FILE_MAP
do
    FILE_ARRAY=($FILE_MAP)
    BASENAME=$(basename "${FILE_ARRAY[0]}" ".txz")
    MakeDb.py imgt -i "${FILE_ARRAY[0]}" -s "${FILE_ARRAY[1]}" \
        --scores --regions --junction --outdir $2/db
    ParseDb.py select -d $2/db/"${BASENAME}_db-pass.tab" -f FUNCTIONAL -u T \
        --outname "${BASENAME}_functional" --outdir $2/db
    ParseDb.py select -d $2/db/"${BASENAME}_functional_parse-select.tab" -f V_CALL J_CALL \
        -u "IGH" --logic all --regex --outname "${BASENAME}_functional-heavy" \
        --outdir $2/functional
    ParseDb.py select -d $2/db/"${BASENAME}_functional_parse-select.tab" -f V_CALL J_CALL \
        -u "IG[LK]" --logic all --regex --outname "${BASENAME}_functional-light" \
        --outdir $2/functional
    rm $2/db/"${BASENAME}_functional_parse-select.tab"
done < $1

