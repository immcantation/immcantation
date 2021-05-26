#!/usr/bin/env python3
"""
Convert IMGT germline fastq files to cellranger-mkvdjref input
"""
# Imports
import argparse
import os
import sys
from pkg_resources import parse_version
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from presto.Annotation import convertIMGTHeader
from changeo.Gene import getGene, getCGene, getAlleleNumber, getLocus

# Check requirements
if parse_version(Bio.__version__) < parse_version('1.71'):
    sys.exit('biopython >= 1.71 is required. Installed version is %s' % Bio.__version__)


def convertHeader(seq, ident=1, pseudogene=False):
    """

    Arguments:
        seq (Bio.Seq): Seq object containing the IMGT formatted header
        ident (int): unique identifier for the record.

    Returns:
        tuple: tuple of strings defining the (gene name, cellranger-mkvdjref header)
    """
    # Parse region_type field
    def _region(x):
        v_regions = ['L-PART1+V-EXON', 'L-PART1+L-PART2+V-REGION', 'V-REGION']
        c_regions = ['CH1', 'CL', 'EX1', 'M1', 'M2', 'C-REGION']
        if x['REGION'] in v_regions:
            return 'L-REGION+V-REGION'
        elif x['REGION'].split('+')[0] in c_regions:
            return 'C-REGION'
        else:
            return x['REGION']

    # Parse gene
    def _gene(x):
        if _region(x) == 'C-REGION':
            return getCGene(x['ID'])
        else:
            return x['ID']

    # Parse isotype field
    def _isotype(x):
        if _region(x) == 'C-REGION' and getLocus(x['ID']) == 'IGH':
            return getCGene(x['ID'])[3:]
        else:
            return 'None'

    # Extract header into dictionary
    imgt = convertIMGTHeader(seq.description)

    # Build cellranger-mkvdjref header
    seq_gene = imgt['ID']
    # print(imgt)
    # print(seq_gene)
    seq_ann = [imgt['ACCESSION'],
               _gene(imgt),
               _region(imgt),
               getLocus(seq_gene)[0:2],
               getLocus(seq_gene),
               _isotype(imgt),
               getAlleleNumber(seq_gene)]

    header = '%i|%s %s' % (ident, seq_gene, '|'.join(seq_ann))
    # print(header)

    # Return valid sequences
    if 'N' in seq.seq.upper():
        # print(seq.seq)
        return seq_gene, None
    if not pseudogene and imgt['FUNCTIONALITY'] not in ('F', 'ORF'):
        # print(imgt)
        return seq_gene, None
    else:
        return seq_gene, header


def main(input_files, output_file, pseudogene=False):
    """
    Convert IMGT germline fastq files to cellranger-mkvdjref input

    Arguments:
        input_files (list):  list of input fasta files.
        output_file (str):  output fasta file.

    Examples:
    >id|display_name record_id|gene_name|region_type|chain_type|chain|isotype|allele_name

    >1|TRAV1*01 AF259072|TRAV1|L-REGION+V-REGION|TR|TRA|None|01
    >1005|IGHG4*04 AL928742|IGHG4|C-REGION|IG|IGH|G4|04

    id              Unique integer ID for this feature.
    display_name    This is used when displaying the segment in, e.g., Loupe V(D)J Browser.
    record_id       Describes the accession ID of the source molecule. Unused.
    gene_name       The name of the V(D)J gene, e.g. TRBV2-1.
    region_type     The only used values are L-REGION+V-REGION, D-REGION, J-REGION, and C-REGION.
    chain_type      Specifies whether this is a T- or B- cell receptor chain. The only used values are TR and IG.
    isotype         Specifies the class of heavy chain constant region; set to None if not applicable.
    allele_name     The identifier for the allele, e.g. 01 for TRBV2-1*01, or None if no allele is to be specified.
    """

    # Iterate over IMGT input files, convert headers, and strip gaps
    name_set = set()
    seq_list = list()
    i = 1
    for f in input_files:
        for rec in Bio.SeqIO.parse(f, 'fasta'):
            # Convert header and validate
            name, header = convertHeader(rec, i)
            if header is None:  continue

            # Append sequence record
            rec = SeqRecord(rec.seq.ungap('.').upper(), id=header, name=header, description='')
            # print(rec)
            if name not in name_set:
                name_set.add(name)
                seq_list.append(rec)
                i += 1

    # Write to mkvdjref input file
    with open(output_file, 'w') as out_handle:
        SeqIO.write(seq_list, out_handle, format='fasta-2line')


if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    
    Examples:
        imgt2cellranger.py -i imgt/rat/leader_vexon/*IG?V.fasta imgt/rat/vdj/*IG?[DJ].fasta imgt/rat/constant/*IG?C.fasta \
            -o Rattus_norvegicus_mkvdjref_input.fasta
        
        imgt2cellranger.py -i imgt/rat/vdj/*IG?[VDJ].fasta imgt/rat/constant/*IG?C.fasta \
            -o Rattus_norvegicus_mkvdjref_input_sans-leader.fasta
    """
    parser = argparse.ArgumentParser(description='Convert IMGT germline fastq files to cellranger-mkvdjref input')
    parser.add_argument('-i', dest='input_files', required=True, nargs='+',
                        help='List of input IMGT reference fastq files.')
    parser.add_argument('-o', dest='output_file', required=True,
                        help='Output fasta file.')
    parser.add_argument('-p', dest='pseudogene', action='store_true',
                        help='If specified retain pseudogenes in the output. By default, only functional and ORF genes that contain no N characters are retained.')
    args = parser.parse_args()

    # Check that files exist
    for f in args.input_files:
        if not os.path.isfile(f):  sys.exit('File %s does not exist.' % f)

    main(**args.__dict__)
