Docker Container
================================================================================

We have provided a complete installation of the Immcantation framework, its
dependencies, accessory scripts, and IgBLAST in a
`Docker <http://www.docker.com>`__. The image also includes both the IgBLAST and
IMGT reference germline sets, as well as several template pipeline scripts.
The image is available on docker hub at:

`kleinstein/immcantation <https://hub.docker.com/r/kleinstein/immcantation/>`__

Images are versioned through tags with images containing official releases
denoted by meta-version numbers (eg, ``1.0.0``). The ``devel`` tag denoted the
latest development (unstabled) builds.

Getting the Container
--------------------------------------------------------------------------------

Requires an installation of Docker 1.9+ or Singularity 2.3+.

Docker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    # Pull release version 1.0.0
    docker pull kleinstein/immcantation:1.0.0

    # Pull the latest development build
    docker pull kleinstein/immcantation:devel


Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    # Pull release version 1.0.0
    IMAGE="immcantation-1.0.0.img"
    singularity create --size 4000 $IMAGE
    singularity import $IMAGE docker://kleinstein/immcantation:1.0.0


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

The following accessory scripts are found in ``/usr/local/bin``:

* ``fastq2fasta.py``:  Simple FASTQ to FASTA conversion.
* ``fetch_igblastdb.sh``:  Downloads the IgBLAST reference database.
* ``fetch_imgtdb.sh``:  Downloads the IMGT reference database.
* ``imgt2igblast.sh``:  Imports the IMGT reference database into IgBLAST.
* ``run_igblast.sh``:  Simple IgBLAST wrapper for running IgBLAST with
  the required arguments using the IMGT reference database.

Template Pipeline Scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following pipeline templates are found in ``/usr/local/bin``:

* ``presto-abseq``:  A start to finish pRESTO processing script for AbSeq data.
* ``changeo-igblast``:  Performs V(D)J alignment using IgBLAST and
  post-processes the output into the Change-O data standard.
* ``changeo-clone``:  Assigns Ig sequences into clonally related lineages and
  builds full germline sequences.
* ``shazam-threshold``:  Performs automated detection of the clonal assignment threshold.
* ``tigger-genotype``:  Infers V segment genotypes using TIgGER.

Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``/usr/local/share/germlines/imgt/IMGT.yaml``:  Information about the downloaded
  IMGT reference sequences.
* ``/usr/local/share/germlines/imgt/<species>/vdj``:  IMGT-gapped V(D)J reference
  sequences.
* ``/usr/local/share/igblast``:  IgBLAST data directory.
* ``/usr/local/share/igblast/fasta``:  Ungapped IMGT references sequences with
  IGH/IGL/IGL and TRA/TRB/TRG/TRD combined into single files, respectively.


Using the Container
--------------------------------------------------------------------------------

Sharing files between the host operating system and the container requires you
to bind one of the container's mount points to a folder on the host using the
``-v`` argument to ``docker`` or the ``-B`` argument to ``singularity``.
There are three available mount points defined in the container::

    /data
    /scratch
    /software

To invoke a shell session inside the container with ``$HOME/project`` mounted to
``/data``::

    # Docker command
    docker run -it -v $HOME/project:data:z kleinstein/immcantation:1.0.0 bash

    # Singularity command
    singularity shell -B $HOME/project:/data immcantation-1.0.0.img

Note, the ``:z`` in the ``-v`` argument of the ``docker`` command is essential.

To execute a specific command::

    # Docker command
    docker run -v $HOME/project:data:z kleinstein/immcantation:1.0.0 versions report

    # Singularity command
    singularity exec -B $HOME/project:/data immcantation-1.0.0.img versions report

In this case, we are executing the ``versions report`` command which will inspect
the installed software versions and print them to standard output.


Running the Template Pipeline Scripts
--------------------------------------------------------------------------------

You can always run your own pipeline scripts through the container, but the
container also includes a set of predefined pipeline scripts that can be run as
is or extended to your needs. Each pipeline script has a ``-h`` argument which
will explain its use. The available pipelines are:

* ``presto-abseq``
* ``changeo-igblast``
* ``changeo-clone``
* ``tigger-genotype``
* ``shazam-threshold``

All template pipeline scripts can be found in ``/usr/local/bin``.

pRESTO pipeline for AbSeq data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A start to finish pRESTO processing script for AbSeq data. Primer sequences are
available from the Immcantation repository under
`protocols/AbSeq <https://bitbucket.org/kleinstein/immcantation/src/tip/protocols/AbSeq>`__

Arguments:
   -1  Read 1 FASTQ sequence file (sequence beginning with the C-region or J-segment).
   -2  Read 2 FASTQ sequence file (sequence beginning with the leader or V-segment).
   -j  Read 1 FASTA primer sequences (C-region or J-segment).
   -v  Read 2 FASTA primer sequences (template switch or V-segment).
   -c  C-region FASTA sequences for the C-region internal to the primer.
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

.. code-block:: shell

    # Arguments
    DATA_DIR=~/project
    READS_R1=/data/raw/sample_R1.fastq
    READS_R2=/data/raw/sample_R2.fastq
    PRIMERS_R1=/data/primers/AbSeqV3_R1_Human_IG_Primers.fasta
    PRIMERS_R2=/data/primers/AbSeqV3_R2_TS.fasta
    CREGION=/data/primers/AbSeqV3_Human_InternalCRegion.fasta
    YAML=/data/sample.yaml
    SAMPLE_NAME=sample
    OUT_DIR=/data/presto/sample
    NPROC=4

    # Docker command
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:1.0.0 presto-abseq \
        -1 $READS_R1 -2 $READS_R2 -j $PRIMERS_R1 -v $PRIMERS_R2 \
        -c $CREGION -y $YAML -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_presto.out

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-1.0.0.img presto-abseq \
        -1 $READS_R1 -2 $READS_R2 -j $PRIMERS_R1 -v $PRIMERS_R2 \
        -c $CREGION -y $YAML -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_presto.out

IgBLAST pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Performs V(D)J alignment using IgBLAST and post-processes the output into the
Change-O data standard.

Arguments:
   -s  FASTA or FASTQ sequence file.
   -r  Directory containing IMGT-gapped reference germlines.
       Defaults to /usr/local/share/germlines/imgt/human/vdj.
   -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file
       and auxillary_data directories. Defaults to /usr/local/share/igblast.
   -n  Sample name or run identifier which will be used as the output file prefix.
       Defaults to a truncated version of the read 1 filename.
   -o  Output directory.
       Defaults to the sample name.
   -p  Number of subprocesses for multiprocessing tools.
       Defaults to the available processing units.
   -h  Display help.

.. code-block:: shell

    # Arguments
    DATA_DIR=~/project
    READS=/data/presto/sample/sample-final_collapse-unique_atleast-2.fastq
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:1.0.0 changeo-igblast \
        -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_igblast.out

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-1.0.0.img changeo-igblast \
        -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_igblast.out

Genotyping pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Infers V segment genotypes using TIgGER.

Arguments:
   -d  Change-O formatted TSV (TAB) file.
   -r  FASTA file containing IMGT-gapped V segment reference germlines.
       Defaults to /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta.
   -n  Sample name or run identifier which will be used as the output file prefix.
       Defaults to a truncated version of the input filename.
   -o  Output directory.
       Defaults to current directory.
   -p  Number of subprocesses for multiprocessing tools.
       Defaults to the available processing units.
   -h  Display help.

.. code-block:: shell

    # Arguments
    DATA_DIR=~/project
    DB=/data/changeo/sample/sample_db-pass.tab
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:1.0.0 tigger-genotype \
        -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_genotype.out

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-1.0.0.img tigger-genotype \
        -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_genotype.out

Clonal threshold inferrence pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Performs automated detection of the clonal assignment threshold.

Arguments:
   -d  Change-O formatted TSV (TAB) file.
   -m  Method.
       Defaults to gmm.
   -n  Sample name or run identifier which will be used as the output file prefix.
       Defaults to a truncated version of the input filename.
   -o  Output directory.
       Defaults to current directory.
   -p  Number of subprocesses for multiprocessing tools.
       Defaults to the available processing units.
   -h  Display help.

.. code-block:: shell

    # Arguments
    DATA_DIR=~/project
    DB=/data/changeo/sample/sample_genotyped.tab
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:1.0.0 shazam-threshold \
        -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_threshold.out

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-1.0.0.img shazam-threshold \
        -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_threshold.out

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
   -h  Display help.

.. code-block:: shell

    # Arguments
    DATA_DIR=~/project
    DB=/data/changeo/sample/sample_genotyped.tab
    DIST=0.15
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image
    docker run -v $DATA_DIR:/data:z kleinstein/immcantation:1.0.0 changeo-clone \
        -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_clone.out

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation-1.0.0.img changeo-clone \
        -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC \
        | tee run_clone.out