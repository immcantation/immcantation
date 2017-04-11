#!/usr/bin/env python3
"""
Converts FASTQ to FASTA
"""
from os import path
from sys import argv
from Bio import SeqIO

infile = argv[1]
outfile = '%s.fasta' % path.splitext(infile)[0]

with open(outfile, 'w') as out_handle:
    seq = SeqIO.parse(infile, 'fastq')
    writer = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)
    writer.write_file(seq)

print(outfile)