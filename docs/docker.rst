Docker Container
================================================================================

A `Docker <http://www.docker.com>`__ image containing a complete installation
of the Immcantation framework, its dependencies, accessory scripts, and IgBLAST
with embedded IMGT reference germlines is available from docker hub under:

`kleinstein/immcantation <https://hub.docker.com/r/kleinstein/immcantation/>`__

Versioned images (tags ``x.y.z``) contained the release builds for the Immcantation
framework. The ``devel`` tag contains the latest development (unstabled) builds.


Getting the Container
--------------------------------------------------------------------------------

Requires and installation of Docker 1.9+ or Singularity 2.3+.

Docker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    docker pull kleinstein/immcantation:devel

Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    IMAGE="immcantation-devel.img"
    singularity create --size 4000 $IMAGE
    singularity import $IMAGE docker://kleinstein/immcantation:devel


What's in the Container
--------------------------------------------------------------------------------

Immcantation Tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `pRESTO <http://presto.readthedocs.io>`__
* `Change-O <http://changeo.readthedocs.io>`__
* `Alakazam <http://alakazam.readthedocs.io>`__
* `SHazaM <http://shazam.readthedocs.io>`__
* `TIgGER <http://tigger.readthedocs.io>`__
* `prestoR <http://bitbucket.org/javh/prototype-prestor>`__

Third Party Tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `muscle <http://www.drive5.com/muscle>`__
* `vsearch <http://github.com/torognes/vsearch>`__
* `BLAST <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`__
* `IgBLAST <https://www.ncbi.nlm.nih.gov/igblast>`__
* `PHYLIP <http://evolution.gs.washington.edu/phylip>`__

Accessory Scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``fastq2fasta.py``:  Simple FASTQ to FASTA conversion.
* ``fetch_igblastdb.sh``:  Downloads the IgBLAST reference database.
* ``fetch_imgtdb.sh``:  Downloads the IMGT reference database.
* ``imgt2igblast.sh``:  Imports the IMGT reference database into IgBLAST.
* ``run_igblast.sh``:  Simple IgBLAST wrapping for running IgBLAST with
  the required arguments and reference database.

Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``/usr/local/share/germlines/imgt/IMGT.yaml``:  Information about the downloaded
  IMGT reference sequences.
* ``/usr/local/share/germlines/imgt/<species>/vdj``:  IMGT-gapped V(D)J reference
  sequences.
* ``/usr/local/share/igblast``:  IgBLAST data directory.
* ``/usr/local/share/igblast/fasta``:  Ungapped IMGT references sequences with
  IGH/IGL/IGL and TRA/TRB/TRG/TRD combined into single files, respectively.

Pipelines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``presto-abseq``:  Start to finish pRESTO processing script for AbSeq data.
* ``changeo-igblast``:  IgBLAST alignment and post processing to the Change-O
  data standard.
* ``changeo-clone``:  Assign Ig sequences into clonal related lineages and
  build full germline sequences.
* ``shazam-threshold``:  Automated detection of clonal assignment threshold.
* ``tigger-genotype``:  V segment genotyping via TIgGER.


Using the Container
--------------------------------------------------------------------------------

Sharing files between the host OS and the container requires you to bind one
of the container's mount points to a folder on the host using the ``-v``
argument to ``docker`` or the ``-B`` argument to ``singularity``.
There are three mount points in the container::

    /data
    /scratch
    /software

To invoke a shell session with ``$HOME/project`` mounted to ``/data``::

    docker run -it -v $HOME/project:data:z kleinstein/immcantation:devel bash
    singularity shell -B $HOME/project:/scratch immcantation-devel.img

Note, the ``:z`` in the ``-v`` argument of the ``docker`` command is essential.

To execute a specific command::

    docker run -v $HOME/project:data:z kleinstein/immcantation:devel versions report
    singularity exec -B $HOME/project:/scratch immcantation-devel.img versions report

In this case, we are executing the ``versions report`` command which will inspect
the installed software versions and print them.

Embedded Pipelines
--------------------------------------------------------------------------------

You can always run your own pipeline scripts through the container, but the
contained also includes a set of predefined pipeline scripts that can be run as
is or extended to your needs. Each pipeline script has a ``-h`` argument which
will explain its use. The available pipelines are:

* ``presto-abseq``
* ``changeo-igblast``
* ``changeo-clone``
* ``tigger-genotype``
* ``shazam-threshold``

Run the pRESTO pipeline for AbSeq data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    # Arguments
    DATA_DIR=~/project
    READS_R1=/data/raw/sample_R1.fastq
    READS_R1=/data/raw/sample_R2.fastq
    PRIMERS_R1=/data/primers/AbSeqV3_Human_R1CPrimers.fasta
    PRIMERS_R2=/data/primers/AbSeqV3_Human_R2TSPrimers.fasta
    CREGION=/data/primers/AbSeqV3_Human_InternalCRegion.fasta
    YAML=/data/sample.yaml
    SAMPLE_NAME=sample
    OUT_DIR=/data/presto/sample
    NPROC=4

    # Run pipeline in docker image.
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:devel presto-abseq \
        -1 $READS_R1 -2 $READS_R2 -j $PRIMERS_R1 -v $PRIMERS_R2 \
        -c $CREGION -y $YAML -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_presto.out

Run the IgBLAST pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    # Arguments
    DATA_DIR=~/project
    READS=/data/presto/sample-final_collapse-unique_atleast-2.fastq
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image.
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:devel changeo-igblast \
        -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_igblast.out

Run the genotyping pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    # Arguments
    DATA_DIR=~/project
    DB=/data/changeo/sample/sample_db-pass.tab
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image.
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:devel tigger-genotype \
        -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_genotype.out

Run the clonal threshold inferrence pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    # Arguments
    DATA_DIR=~/project
    DB=/data/changeo/sample/sample_db-pass.tab
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image.
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:devel shazam-threshold \
        -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_threshold.out

Run the clonal assignment pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    # Arguments
    DATA_DIR=~/project
    DB=/data/changeo/sample/sample_db-pass.tab
    DIST=0.15
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image.
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:devel changeo-clone \
        -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_clone.out
