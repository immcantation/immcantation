Docker Container Installation and Overview
================================================================================

We have provided a complete installation of the Immcantation framework, its
dependencies, accessory scripts, and IgBLAST in a
`Docker <http://www.docker.com>`__. The image also includes both the IgBLAST and
IMGT reference germline sets, as well as several template pipeline scripts.
The image is available on docker hub at:

`kleinstein/immcantation <https://hub.docker.com/r/kleinstein/immcantation/>`__

Images are versioned through tags with images containing official releases
denoted by meta-version numbers (``x.y.z``). The ``devel`` tag denotes the
latest development (unstabled) builds.

Getting the Container
--------------------------------------------------------------------------------

Requires an installation of Docker 1.9+ or Singularity 2.3+.

Docker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. parsed-literal::

    # Pull release version |docker-version|
    docker pull kleinstein/immcantation:|docker-version|

    # Pull the latest development build
    docker pull kleinstein/immcantation:devel


Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. parsed-literal::

    # Pull release version |docker-version|
    IMAGE="immcantation-|docker-version|.img"
    singularity create --size 6000 $IMAGE
    singularity import $IMAGE docker://kleinstein/immcantation:|docker-version|


What's in the Container
--------------------------------------------------------------------------------

Immcantation Tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `pRESTO <http://presto.readthedocs.io>`__
* `Change-O <http://changeo.readthedocs.io>`__
* `Alakazam <http://alakazam.readthedocs.io>`__
* `SHazaM <http://shazam.readthedocs.io>`__
* `TIgGER <http://tigger.readthedocs.io>`__
* `RDI <http://rdi.readthedocs.io>`__
* `prestoR <http://bitbucket.org/javh/prototype-prestor>`__
* `SCOPe <http://bitbucket.org/kleinstein/scope>`__

Third Party Tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `muscle <http://www.drive5.com/muscle>`__
* `vsearch <http://github.com/torognes/vsearch>`__
* `BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`__
* `IgBLAST <https://www.ncbi.nlm.nih.gov/igblast>`__
* `IgPhyML <https://bitbucket.org/kbhoehn/igphyml>`__
* `PHYLIP <http://evolution.gs.washington.edu/phylip>`__

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
run_igblast.sh
    Simple IgBLAST wrapper for running IgBLAST with the required arguments
    using the IMGT reference database.

Template Pipeline Scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following pipeline templates are found in ``/usr/local/bin``:

presto-abseq
    A start to finish pRESTO processing script for AbSeq data.
changeo-igblast
    Performs V(D)J alignment using IgBLAST and post-processes the output into
    the Change-O data standard.
changeo-clone
    Assigns Ig sequences into clonally related lineages and builds full
    germline sequences.
shazam-threshold
    Performs automated detection of the clonal assignment threshold.
tigger-genotype
    Infers V segment genotypes using TIgGER.
preprocess-phix
    Removes PhiX reads from raw data files.

Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``/usr/local/share/germlines/imgt/IMGT.yaml``
    Information about the downloaded IMGT reference sequences.
``/usr/local/share/germlines/imgt/<species>/vdj/``
    Directory containing IMGT-gapped V(D)J reference sequences in FASTA format.
``/usr/local/share/igblast/``
    IgBLAST data directory.
``/usr/local/share/igblast/fasta/``
    Directory containing ungapped IMGT references sequences with IGH/IGL/IGL and
    TRA/TRB/TRG/TRD combined into single FASTA files, respectively.
``/usr/local/share/protocols``
    Directory containing primer, template switch and internal constant region
    sequences for various experimental protocols in FASTA format.
