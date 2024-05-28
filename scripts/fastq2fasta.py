#!/usr/bin/env python3
"""
Converts FASTQ to FASTA
"""
import Bio
from packaging.version import Version
from os import path
from sys import argv
from Bio import SeqIO

in_file = argv[1]
out_file = path.split(in_file)[1]
out_file = '%s.fasta' % path.splitext(out_file)[0]

with open(out_file, 'w') as out_handle:
    records = SeqIO.parse(in_file, 'fastq')
    if Version(Bio.__version__) >= Version('1.71'):
        # Biopython >= v1.71
        SeqIO.write(records, out_handle, format='fasta-2line')
    else:
        # Biopython < v1.71
        writer = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)
        writer.write_file(records)

print(out_file)
