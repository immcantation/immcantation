What's in the /suite Container
--------------------------------------------------------------------------------

Immcantation Tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+ `pRESTO <https://presto.readthedocs.io>`__
+ `Change-O <https://changeo.readthedocs.io>`__
+ `Alakazam <https://alakazam.readthedocs.io>`__
+ `SHazaM <https://shazam.readthedocs.io>`__
+ `TIgGER <https://tigger.readthedocs.io>`__
+ `SCOPer <https://scoper.readthedocs.io>`__
+ `dowser <https://dowser.readthedocs.io>`__
+ `RDI <https://rdi.readthedocs.io>`__
+ `prestoR <https://bitbucket.org/kleinstein/prestor>`__

Third Party Tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+ `muscle <http://www.drive5.com/muscle>`__
+ `vsearch <http://github.com/torognes/vsearch>`__
+ `CD-HIT <http://weizhongli-lab.org/cd-hit>`__
+ `BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`__
+ `IgBLAST <https://www.ncbi.nlm.nih.gov/igblast>`__
+ `IgPhyML <https://bitbucket.org/kleinstein/igphyml>`__
+ `PHYLIP <http://evolution.gs.washington.edu/phylip>`__
+ `tbl2asn <https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2>`__

Example Pipeline Scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+ :ref:`PhiXPipeline`
+ :ref:`AbSeqPipeline`
+ :ref:`ClontechPipeline`
+ :ref:`ClontechPipelineUMI`
+ :ref:`10XPipeline`
+ :ref:`IgBLASTPipeline`
+ :ref:`ClonePipeline`
+ :ref:`ThresholdPipeline`
+ :ref:`GenotypePipeline`

Accessory Scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following accessory scripts are found in ``/usr/local/bin``:

fastq2fasta.py
    Simple FASTQ to FASTA conversion.
fetch_phix.sh
    Downloads the PhiX174 reference genome.
fetch_igblastdb.sh
    Downloads the IgBLAST reference database.
fetch_imgtdb.sh
    Downloads the IMGT reference database.
imgt2igblast.sh
    Imports the IMGT reference database into IgBLAST.
imgt2cellranger.py
    Converts the IMGT fasta germline reference files to the input required by
    cellranger-mkvdjref.

Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``/usr/local/share/germlines/imgt/IMGT.yaml``
    Information about the downloaded IMGT reference sequences.
``/usr/local/share/germlines/imgt/<species>/vdj``
    Directory containing IMGT-gapped V(D)J reference sequences in FASTA format.
``/usr/local/share/igblast``
    IgBLAST data directory.
``/usr/local/share/igblast/fasta``
    Directory containing ungapped IMGT references sequences with IGH/IGK/IGL and
    TRA/TRB/TRG/TRD combined into single FASTA files, respectively.
``/usr/local/share/protocols``
    Directory containing primer, template switch and internal constant region
    sequences for various experimental protocols in FASTA format.
