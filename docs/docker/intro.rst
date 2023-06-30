.. _DockerIntro:

Docker Container Installation and Overview
================================================================================

We have provided a complete installation of the Immcantation framework, its
dependencies, accessory scripts, and IgBLAST in a
`Docker container <http://www.docker.com>`__. The `container <https://hub.docker.com/r/immcantation/suite/>`__ also includes both the IgBLAST and
IMGT reference germline sets, as well as several template pipeline scripts.

We currently have four containers available on `DockerHub <https://hub.docker.com/r/immcantation/>`__:

+---------------------------------------+-----------------------------------------------------------------------------------------+
| Name                                  | Contents                                                                                |
+=======================================+=========================================================================================+
| immcantation/suite                    | Immcantation suite, supporting applications and databases.                              |
+---------------------------------------+-----------------------------------------------------------------------------------------+
| immcantation/lab                      | Immcantation tutorial materials. Only for training, not to be used in production.       |
+---------------------------------------+-----------------------------------------------------------------------------------------+
| immcantation/base                     | Base image for Immcantation development builds.                                         |
+---------------------------------------+-----------------------------------------------------------------------------------------+
| immcantation/test                     | Immcantation unit test image.                                                           |
+---------------------------------------+-----------------------------------------------------------------------------------------+

**For tutorial purposes**, use immcantation/lab (be sure to replace ``suite`` with ``lab`` in the following code chunks) and follow the directions :ref:`here <DockerGuideTutorials>`. For all Immcantation uses, use immcantation/suite.
Note that containers are versioned through tags with containers containing official releases
denoted by meta-version numbers (``x.y.z``). The ``devel`` tag denotes the
latest development (unstable) builds.

Getting the /suite Container
--------------------------------------------------------------------------------

Requires an installation of Docker 1.9+ or Singularity 2.3+.

Docker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. parsed-literal::

    # Pull release version |docker-version|
    docker pull immcantation/suite:|docker-version|

    # Pull the latest development build
    docker pull immcantation/suite:devel


Our containers are Linux-based, so if you are using a Windows computer,
please make sure that you are using Linux containers and not Windows containers
(this can be changed in `Docker Desktop <https://www.docker.com/products/docker-desktop/>`_
and won't affect your existing containers).


Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. parsed-literal::

    # Pull release version |docker-version|
    IMAGE="immcantation_suite-|docker-version|.sif"
    singularity build $IMAGE docker://immcantation/suite:|docker-version|

The instructions to use containers from `Docker Hub <https://hub.docker.com/>`_
with Singularity can be slightly different for different versions of Singularity.
If the command shown above doesn't work for you, please visit
`Singularity Documentation <https://www.sylabs.io/docs/>`_ and look for the
specific command for your Singularity version under *Build a container*.


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

Template Pipeline Scripts
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
