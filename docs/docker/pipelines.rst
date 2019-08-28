.. _PipelineScripts:

Pipeline Templates
================================================================================

You can always run your own pipeline scripts through the container, but the
container also includes a set of predefined pipeline scripts that can be run as
is or extended to your needs. Each pipeline script has a ``-h`` argument which
will explain its use. The available pipelines are:

* preprocess-phix
* presto-abseq
* presto-clontech
* changeo-10x
* changeo-igblast
* tigger-genotype
* shazam-threshold
* changeo-clone

All template pipeline scripts can be found in ``/usr/local/bin``.

PhiX cleaning pipeline
--------------------------------------------------------------------------------

Removes reads from a sequence file that align against the PhiX174 reference
genome.

.. include:: ../_include/usage.rst
    :start-after: Start preprocess-phix
    :end-before: End preprocess-phix

**Example: preprocess-phix**

.. parsed-literal::

    # Arguments
    DATA_DIR=~/project
    READS=/data/raw/sample.fastq
    OUT_DIR=/data/presto/sample
    NPROC=4

    # Run pipeline in docker image
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:|docker-version| \\
        preprocess-phix -s $READS -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.sif \\
        preprocess-phix -s $READS -o $OUT_DIR -p $NPROC

.. note::

    The PhiX cleaning pipeline will convert the sequence headers to
    the pRESTO format. Thus, if the ``nophix`` output file is provided as
    input to the ``presto-abseq`` pipeline script you must pass the argument
    ``-x presto`` to ``presto-abseq``, which will tell the
    script that the input headers are in pRESTO format (rather than the
    Illumina format).


NEB AbSeq protocol pRESTO pipeline
--------------------------------------------------------------------------------

A start to finish pRESTO processing script for AbSeq data. Primer sequences are
available from the Immcantation repository under
`protocols/AbSeq <https://bitbucket.org/kleinstein/immcantation/src/tip/protocols/AbSeq>`__
or inside the container under ``/usr/local/share/protocols/AbSeq``.

.. include:: ../_include/usage.rst
    :start-after: Start presto-abseq
    :end-before: End presto-abseq

One of the requirements for generating the report at the end of the pRESTO pipeline is a YAML
file containing information about the data and processing. Valid fields are shown in the example
``sample.yaml`` below, although no fields are strictly required:

**sample.yaml**

.. parsed-literal::

    title: "pRESTO Report: CD27+ B cells from subject HD1"
    author: "Your Name"
    version: "0.5.4"
    description: "Memory B cells (CD27+)."
    sample: "HD1"
    run: "ABC123"
    date: "Today"

**Example: presto-abseq**

.. parsed-literal::

    # Arguments
    DATA_DIR=~/project
    READS_R1=/data/raw/sample_R1.fastq
    READS_R2=/data/raw/sample_R2.fastq
    YAML=/data/sample.yaml
    SAMPLE_NAME=sample
    OUT_DIR=/data/presto/sample
    NPROC=4

    # Docker command
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:|docker-version| \\
        presto-abseq -1 $READS_R1 -2 $READS_R2 -y $YAML \\
        -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.sif \\
        presto-abseq -1 $READS_R1 -2 $READS_R2 -y $YAML \\
        -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC


Takara Bio / Clontech SMARTer protocol pRESTO pipeline
--------------------------------------------------------------------------------

A start to finish pRESTO processing script for Takara Bio / Clontech SMARTer kit
data. C-regions are assigned using the universal C-region primer sequences are
available from the Immcantation repository under
`protocols/Universal <https://bitbucket.org/kleinstein/immcantation/src/tip/protocols/Universal>`__
or inside the container under ``/usr/local/share/protocols/Universal``.

.. include:: ../_include/usage.rst
    :start-after: Start presto-clontech
    :end-before: End presto-clontech

**Example: presto-clontech**

.. parsed-literal::

    # Arguments
    DATA_DIR=~/project
    READS_R1=/data/raw/sample_R1.fastq
    READS_R2=/data/raw/sample_R2.fastq
    CREGION=/usr/local/share/protocols/Universal/Human_IG_CRegion_RC.fasta
    VREF=/usr/local/share/igblast/fasta/imgt_human_ig_v.fasta
    SAMPLE_NAME=sample
    OUT_DIR=/data/presto/sample
    NPROC=4

    # Docker command
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:|docker-version| \\
        presto-clontech -1 $READS_R1 -2 $READS_R2 -j $CREGION -r $VREF \\
        -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.sif \\
        presto-abseq -1 $READS_R1 -2 $READS_R2 -j $CREGION -r $VREF \\
        -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC


10X Genomics V(D)J annotation pipeline
--------------------------------------------------------------------------------

Assigns new annotations and infers clonal relationships to 10X Genomics
single-cell V(D)J data output by Cell Ranger.

.. include:: ../_include/usage.rst
    :start-after: Start changeo-10x
    :end-before: End changeo-10x

**Example: changeo-10x**

.. parsed-literal::

    # Arguments
    DATA_DIR=~/project
    READS=/data/raw/sample_filtered_contig.fasta
    ANNOTATIONS=/data/raw/sample_filtered_contig_annotations.csv
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    DIST=auto
    NPROC=4

    # Run pipeline in docker image
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:|docker-version| \\
        changeo-10x -s $READS -a $ANNOTATIONS -x $DIST -n $SAMPLE_NAME \\
        -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.sif \\
        changeo-10x -s $READS -a $ANNOTATIONS -x $DIST -n $SAMPLE_NAME \\
        -o $OUT_DIR -p $NPROC


IgBLAST pipeline
--------------------------------------------------------------------------------

Performs V(D)J alignment using IgBLAST and post-processes the output into the
Change-O data standard.

.. include:: ../_include/usage.rst
    :start-after: Start changeo-igblast
    :end-before: End changeo-igblast

**Example: changeo-igblast**

.. parsed-literal::

    # Arguments
    DATA_DIR=~/project
    READS=/data/presto/sample/sample-final_collapse-unique_atleast-2.fastq
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:|docker-version| \\
        changeo-igblast -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.sif \\
        changeo-igblast -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC


Genotyping pipeline
--------------------------------------------------------------------------------

Infers V segment genotypes using TIgGER.

.. include:: ../_include/usage.rst
    :start-after: Start tigger-genotype
    :end-before: End tigger-genotype

**Example: tigger-genotype**

.. parsed-literal::

    # Arguments
    DATA_DIR=~/project
    DB=/data/changeo/sample/sample_db-pass.tab
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:|docker-version| \\
        tigger-genotype -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.sif \\
        tigger-genotype -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

Clonal threshold inference pipeline
--------------------------------------------------------------------------------

Performs automated detection of the clonal assignment threshold.

.. include:: ../_include/usage.rst
    :start-after: Start shazam-threshold
    :end-before: End shazam-threshold

**Example: shazam-threshold**

.. parsed-literal::

    # Arguments
    DATA_DIR=~/project
    DB=/data/changeo/sample/sample_genotyped.tab
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:|docker-version| \\
        shazam-threshold -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.sif \\
        shazam-threshold -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC


Clonal assignment pipeline
--------------------------------------------------------------------------------

Assigns Ig sequences into clonally related lineages and builds full germline
sequences.

.. include:: ../_include/usage.rst
    :start-after: Start changeo-clone
    :end-before: End changeo-clone

**Example: changeo-clone**

.. parsed-literal::

    # Arguments
    DATA_DIR=~/project
    DB=/data/changeo/sample/sample_genotyped.tab
    DIST=0.15
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:|docker-version| \\
        changeo-clone -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.sif \\
        changeo-clone -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC
