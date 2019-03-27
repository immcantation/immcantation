#!/usr/bin/env python3
"""
Bypass MaskPrimers.py align for identification of internal CREGION, adds a CREGION column.
Useful for bypassing MaskPrimers.py in presto by running changeo instead. 

To use: find_internalcregion.py input.tab output.tab AbSeqV3_Human_InternalCRegion.fasta
"""

import sys
import pandas as pd

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from presto.Sequence import localAlignment, maskSeq, getDNAScoreDict
from presto.IO import readPrimerFile
from presto.Defaults import default_primer_max_error, default_primer_max_len, default_barcode_field, \
default_primer_field, default_primer_gap_penalty, default_delimiter
from presto.Annotation import parseAnnotation

input_file = sys.argv[1]
output_file = sys.argv[2]
primer_file = sys.argv[3]

primers = readPrimerFile(primer_file)

def maskInternalCregion(sequence, primers=primers):
    
    # set default parameters (used in MaskPrimers.py)
    max_error=default_primer_max_error
    max_len=default_primer_max_len
    rev_primer=False
    skip_rc=False
    mode='tag'
    barcode=False
    barcode_field=default_barcode_field
    primer_field=default_primer_field
    delimiter=default_delimiter
    gap_penalty=default_primer_gap_penalty
    score_dict=getDNAScoreDict(mask_score=(0, 1), gap_score=(0, 0))
    primers_regex=None
    
    # create a dummy input SeqRecord
    input_seqrecord = SeqRecord(Seq(sequence))

    # run maskprimers on input SeqRecord
    align = localAlignment(input_seqrecord, primers, primers_regex=primers_regex, max_error=max_error,
                       max_len=max_len, rev_primer=rev_primer, skip_rc=skip_rc,
                       gap_penalty=gap_penalty, score_dict=score_dict)
    seqrecord = maskSeq(align, mode=mode, barcode=barcode, barcode_field=barcode_field,
                      primer_field=primer_field, delimiter=delimiter)

    # parse output
    call = parseAnnotation(seqrecord.id)['PRIMER']
    
    return call

# input DF, run MaskPrimers, output to changeo output
input_df = pd.read_csv(input_file, sep = '\t', dtype = 'object')
input_df['CREGION'] = input_df.apply(lambda row: maskInternalCregion(row['SEQUENCE_INPUT'], axis = 1))
input_df.to_csv(output_file , sep = '\t')