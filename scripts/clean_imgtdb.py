#!/usr/bin/env python3
"""
Clean IMGT germline fasta files for IgBLAST database build
"""
import Bio
from packaging.version import Version
from sys import argv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Get input and output file names
in_file = argv[1]
out_file = argv[2]

# Load sequences into memory and process them
name_set = set()
seq_list = list()
seq_list_desc = list()
for rec in SeqIO.parse(in_file, 'fasta'):
    name = rec.description.split('|')[1]
    if name not in name_set:
        name_set.add(name)
        seq = SeqRecord(rec.seq.replace('.', '').upper(), id=name, name=name, description=name)
        seq_desc = SeqRecord(rec.seq.replace('.', '').upper(), id=name, name=name, description=rec.description)
        seq_list.append(seq)
        seq_list_desc.append(seq_desc)
    else:
        # Report duplicated entries (by name) found. Keep the first one and skip the rest.
        # Issue found with constant region sequences, with the same name but different 
        # terminal domains (see descriptions: CH-S vs M).
        # Not relevant for current usage of immcantation, because we don't use sequencing data 
        # from that region.
        # Example:
        # Duplicate sequence name found: IGHA1*01
        # Description: J00220+M60193|IGHA1*01|Homo_sapiens|F|CH1+H+CH2+CH3+M|J00220:g,142..447+662..1021+1244..1575;M60193:g,226..437|1212 nt|1|+1| | | |1212+0=1212| | |
        # Duplicate of sequence with name IGHA1*01 and description: J00220|IGHA1*01|Homo_sapiens|F|CH1+H+CH2+CH3+CH-S|g,142..447+662..1021+1244..1635|1059 nt|1|+1| | | |1059+0=1059| | |
        # Skipping duplicate sequence.
        print(f"Duplicate sequence name found: {name}")
        print(f"Description: {rec.description}")
        print(f"Duplicate of sequence with name {name} and description: {[seq.description for seq in seq_list_desc if seq.name == name][0]}")
        print("Skipping duplicate sequence.\n")

# Overwrite file
with open(out_file, 'w') as out_handle:
    if Version(Bio.__version__) >= Version('1.71'):
        # Biopython >= v1.71
        SeqIO.write(seq_list, out_handle, format='fasta-2line')
    else:
        # Biopython < v1.71
        writer = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)
        writer.write_file(seq_list)
