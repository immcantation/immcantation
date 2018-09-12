Using the Docker Container
================================================================================

Invoking a shell inside the container
--------------------------------------------------------------------------------

To invoke a shell session inside the container:

.. parsed-literal::

    # Docker command
    docker run -it kleinstein/immcantation:|docker-version| bash

    # Singularity command
    singularity shell immcantation-|docker-version|.img

Sharing files with the container
--------------------------------------------------------------------------------

Sharing files between the host operating system and the container requires you
to bind one of the container's mount points to a folder on the host using the
``-v`` argument to ``docker`` or the ``-B`` argument to ``singularity``.
There are four available mount points defined in the container::

    /data
    /scratch
    /software
    /oasis

To invoke a shell session inside the container with ``$HOME/project`` mounted to
``/data``:

.. parsed-literal::

    # Docker command
    docker run -it -v $HOME/project:/data:z kleinstein/immcantation:|docker-version| bash

    # Singularity command
    singularity shell -B $HOME/project:/data immcantation-|docker-version|.img

Note, the ``:z`` in the ``-v`` argument of the ``docker`` command is essential.


Executing a specific command
--------------------------------------------------------------------------------

To execute a specific command inside the container with ``$HOME/project`` mounted to
``/data``:

.. parsed-literal::

    # Docker command
    docker run -v $HOME/project:/data:z kleinstein/immcantation:|docker-version| versions report

    # Singularity command
    singularity exec -B $HOME/project:/data immcantation-|docker-version|.img versions report

In this case, we are executing the ``versions report`` command which will inspect
the installed software versions and print them to standard output.

There is an analagous ``builds report`` command to display the build date and
changesets used during the image build. This is particularly relevant if you
are using the ``kleinstein/immcantation:devel`` development builds.


Running the Template Pipeline Scripts
--------------------------------------------------------------------------------

You can always run your own pipeline scripts through the container, but the
container also includes a set of predefined pipeline scripts that can be run as
is or extended to your needs. Each pipeline script has a ``-h`` argument which
will explain its use. The available pipelines are:

* presto-abseq
* changeo-igblast
* changeo-clone
* tigger-genotype
* shazam-threshold
* preprocess-phix

All template pipeline scripts can be found in ``/usr/local/bin``.

pRESTO pipeline for preprocessing AbSeq data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A start to finish pRESTO processing script for AbSeq data. Primer sequences are
available from the Immcantation repository under
`protocols/AbSeq <https://bitbucket.org/kleinstein/immcantation/src/tip/protocols/AbSeq>`__
or inside the container under ``/usr/local/share/protocols/AbSeq``.

Arguments:
   -1  Read 1 FASTQ sequence file (sequence beginning with the C-region or J-segment).
   -2  Read 2 FASTQ sequence file (sequence beginning with the leader or V-segment).
   -j  Read 1 FASTA primer sequences (C-region or J-segment).
       Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R1_Human_IG_Primers.fasta
   -v  Read 2 FASTA primer sequences (template switch or V-segment).
       Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R2_TS.fasta.
   -c  C-region FASTA sequences for the C-region internal to the primer.
       If unspecified internal C-region alignment is not performed.
   -r  V-segment reference file.
       Defaults to /usr/local/share/germlines/igblast/fasta/imgt_human_ig_v.fasta
   -y  YAML file providing description fields for report generation.
   -n  Sample name or run identifier which will be used as the output file prefix.
       Defaults to a truncated version of the read 1 filename.
   -o  Output directory.
       Defaults to the sample name.
   -p  Number of subprocesses for multiprocessing tools.
       Defaults to the available processing units.
   -h  Display help.

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

**AbSeq preprocessing example**

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
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.img \\
        presto-abseq -1 $READS_R1 -2 $READS_R2 -y $YAML \\
        -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

IgBLAST pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Performs V(D)J alignment using IgBLAST and post-processes the output into the
Change-O data standard.

Arguments:
   -s  FASTA or FASTQ sequence file.
   -r  Directory containing IMGT-gapped reference germlines.
       Defaults to /usr/local/share/germlines/imgt/human/vdj.
   -g  Species name. One of human or mouse. Defaults to human.
   -t  Receptor type. One of ig or tr. Defaults to ig.
   -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file
       and auxillary_data directories. Defaults to /usr/local/share/igblast.
   -n  Sample name or run identifier which will be used as the output file prefix.
       Defaults to a truncated version of the read 1 filename.
   -o  Output directory.
       Defaults to the sample name.
   -p  Number of subprocesses for multiprocessing tools.
       Defaults to the available processing units.
   -f  Specify to filter the output to only productive/functional sequences.
   -h  Display help.

**IgBLAST example**

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
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.img \\
        changeo-igblast -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

Genotyping pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Infers V segment genotypes using TIgGER.

Arguments:
   -d  Change-O formatted TSV (TAB) file.
   -r  FASTA file containing IMGT-gapped V segment reference germlines.
       Defaults to /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta.
   -v  Name of the output field containing genotyped V assignments.
       Defaults to V_CALL_GENOTYPED.
   -n  Sample name or run identifier which will be used as the output file prefix.
       Defaults to a truncated version of the input filename.
   -o  Output directory.
       Defaults to current directory.
   -p  Number of subprocesses for multiprocessing tools.
       Defaults to the available processing units.
   -h  Display help.

**Genotyping example**

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
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.img \\
        tigger-genotype -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

Clonal threshold inferrence pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Performs automated detection of the clonal assignment threshold.

Arguments:
   -d           Change-O formatted TSV (TAB) file.
   -m           Method.
                Defaults to density.
   -n           Sample name or run identifier which will be used as the output file prefix.
                Defaults to a truncated version of the input filename.
   -o           Output directory.
                Defaults to current directory.
   -p           Number of subprocesses for multiprocessing tools.
                Defaults to the available processing units.
   --model      Model when "-m gmm" is specified.
                Defaults to "gamma-gamma".
   --subsample  Number of distances to downsample the data to before threshold calculation.
   --repeats    Number of times to repeat the threshold calculation (with plotting).
   -h           Display help.

**Clonal threshold inferrence example**

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
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.img \\
        shazam-threshold -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

Clonal assignment pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assigns Ig sequences into clonally related lineages and builds full germline
sequences.

Arguments:
   -d  Change-O formatted TSV (TAB) file.
   -x  Distance threshold for clonal assignment.
   -r  Directory containing IMGT-gapped reference germlines.
       Defaults to /usr/local/share/germlines/imgt/human/vdj.
   -n  Sample name or run identifier which will be used as the output file prefix.
       Defaults to a truncated version of the input filename.
   -o  Output directory.
       Defaults to the sample name.
   -p  Number of subprocesses for multiprocessing tools.
       Defaults to the available processing units.
   -a  Specify to clone the full data set.
       By default the data will be filtering to only productive/functional sequences.
   -h  Display help.

**Clonal assignment example**

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
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.img \\
        changeo-clone -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

PhiX cleaning pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Removes reads from a sequence file that align against the PhiX174 reference
genome.

Arguments:
   -s  FASTQ sequence file.
   -r  Directory containing phiX174 reference db.
   -o  Output directory.
       Defaults to the FASTQ file directory.
   -n  Name to use as the output file suffix.
       Defaults to '_nophix'.
   -p  Number of subprocesses for multiprocessing tools.
       Defaults to the available processing units.
   -h  Display help

**PhiX cleaning example**

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
    singularity exec -B $DATA_DIR:/data immcantation-|docker-version|.img \\
        preprocess-phix -s $READS -o $OUT_DIR -p $NPROC
