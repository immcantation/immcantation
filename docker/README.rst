Immcantation Docker Container
=============================

Download
--------

::

    docker pull kleinstein/immcantation:devel

Singularity
-----------

Create a Singularity container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a container from a specific version of Immcantation (do not use
``latest``, this requires ``singularity`` > 2.3, you **do not need
root**.

::

    singularity pull docker://kleinstein/immcantation:3.0.0

Configure bind points
~~~~~~~~~~~~~~~~~~~~~

The container has ``/data`` and ``/scratch`` mountpoints that other used
to make other filesystems (``/home`` is automatically mounted) inside
the container.

Add this variable to your ``.bashrc`` to configure this permanently on
your machine:

::

    export SINGULARITY_BINDPATH="/full/path/to/data/folder:/data,/full/path/to/your/scratch:/scratch"

Interactive shell inside the container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    IMAGE=immcantation-3.0.0.sif
    singularity shell $IMAGE

Launch Immcantation commands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For example:

::

    singularity exec $IMAGE presto-abseq -h

List available pipelines
------------------------

::

    docker run -it kleinstein/immcantation:devel

Specific pipeline help
----------------------

::

    docker run -it kleinstein/immcantation:devel presto-abseq -h

Invoke shell session within image
---------------------------------

::

    DATA_DIR=~/workspace/igpipeline/docker/tmp
    docker run -it -v $DATA_DIR:/data:z kleinstein/immcantation:devel bash

Run pRESTO pipeline for AbSeq v3 data
-------------------------------------

::

    # Arguments
    DATA_DIR=~/workspace/igpipeline/docker/tmp
    READS_R1=/data/raw/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq
    READS_R2=/data/raw/AAYHL_HD13M/MG2v3_HD13M_BC13_AGTCAA_L001_R2_001.fastq
    PRIMERS_R1=/data/primers/AbSeqV3_Human_R1CPrimers.fasta
    PRIMERS_R2=/data/primers/AbSeqV3_Human_R2TSPrimers.fasta
    CREGION=/data/primers/AbSeqV3_Human_InternalCRegion.fasta
    YAML=/data/test.yaml
    SAMPLE_NAME=AAYHL_HD13M
    OUT_DIR=/data/presto/AAYHL_HD13M
    NPROC=4

    # Run pipeline in docker image.
    # Note: mounting a host directory as a volume with the 'z' option is essential.
    docker run -it -v $DATA_DIR:/data:z kleinstein/immcantation:devel presto-abseq \
        -1 $READS_R1 -2 $READS_R2 -j $PRIMERS_R1 -v $PRIMERS_R2 \
        -c $CREGION -y $YAML -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_presto.out

Run IgBLAST pipeline
--------------------

::

    # Arguments
    DATA_DIR=~/workspace/igpipeline/docker/tmp
    READS=/data/presto/AAYHL_HD13M/AAYHL_HD13M-final_collapse-unique_atleast-2.fastq
    SAMPLE_NAME=AAYHL_HD13M
    OUT_DIR=/data/changeo/AAYHL_HD13M
    NPROC=4

    # Run pipeline in docker image.
    # Note: mounting a host directory as a volume with the 'z' option is essential.
    docker run -it -v $DATA_DIR:/data:z kleinstein/immcantation:devel changeo-igblast \
        -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_igblast.out

Run genotyping pipeline
-----------------------

::

    # Arguments
    DATA_DIR=~/workspace/igpipeline/docker/tmp
    DB=/data/changeo/AAYHL_HD13M/AAYHL_HD13M_db-pass.tab
    DIST=0.15
    SAMPLE_NAME=AAYHL_HD13M
    OUT_DIR=/data/changeo/AAYHL_HD13M
    NPROC=4

    # Run pipeline in docker image.
    # Note: mounting a host directory as a volume with the 'z' option is essential.
    docker run -it -v $DATA_DIR:/data:z kleinstein/immcantation:devel tigger-genotype \
        -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_clone.out

Run clonal threshold inferrence
-------------------------------

::

    # Arguments
    DATA_DIR=~/workspace/igpipeline/docker/tmp
    DB=/data/changeo/AAYHL_HD13M/AAYHL_HD13M_db-pass.tab
    SAMPLE_NAME=AAYHL_HD13M
    OUT_DIR=/data/changeo/AAYHL_HD13M
    NPROC=4

    # Run pipeline in docker image.
    # Note: mounting a host directory as a volume with the 'z' option is essential.
    docker run -it -v $DATA_DIR:/data:z kleinstein/immcantation:devel shazam-threshold \
        -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_threshold.out

Run clonal assignment pipeline
------------------------------

::

    # Arguments
    DATA_DIR=~/workspace/igpipeline/docker/tmp
    DB=/data/changeo/AAYHL_HD13M/AAYHL_HD13M_db-pass.tab
    DIST=0.15
    SAMPLE_NAME=AAYHL_HD13M
    OUT_DIR=/data/changeo/AAYHL_HD13M
    NPROC=4

    # Run pipeline in docker image.
    # Note: mounting a host directory as a volume with the 'z' option is essential.
    docker run -it -v $DATA_DIR:/data:z kleinstein/immcantation:devel changeo-clone \
        -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_clone.out
