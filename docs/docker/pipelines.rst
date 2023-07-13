.. _PipelineScripts:

Example pipelines
--------------------------------------------------------------------------------

You can always run your own pipeline scripts through the container, but the
container also includes a set of predefined pipeline scripts that can be run as
is or extended to your needs. Each pipeline script has a ``-h`` argument which
will explain its use. The available pipelines are:

* preprocess-phix
* presto-abseq
* presto-clontech
* presto-clontech-umi
* changeo-10x
* changeo-igblast
* tigger-genotype
* shazam-threshold
* changeo-clone

All example pipeline scripts can be found in ``/usr/local/bin``.


.. _PhiXPipeline:

PhiX cleaning pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    docker run -v $DATA_DIR:/data\:z immcantation/suite:|docker-version| \\
        preprocess-phix -s $READS -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation_suite-|docker-version|.sif \\
        preprocess-phix -s $READS -o $OUT_DIR -p $NPROC

.. note::

    The PhiX cleaning pipeline will convert the sequence headers to
    the pRESTO format. Thus, if the ``nophix`` output file is provided as
    input to the ``presto-abseq`` pipeline script you must pass the argument
    ``-x presto`` to ``presto-abseq``, which will tell the
    script that the input headers are in pRESTO format (rather than the
    Illumina format).


.. _AbSeqPipeline:

NEBNext / AbSeq immune sequencing kit preprocessing pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A start to finish pRESTO processing script for NEBNext / AbSeq immune sequencing data.
An example for human BCR processing is shown below. Primer sequences are available from the
Immcantation repository under `protocols/AbSeq <https://bitbucket.org/kleinstein/immcantation/src/master/protocols/AbSeq>`__
or inside the container under ``/usr/local/share/protocols/AbSeq``. Mouse primers are not supplied.
TCR V gene references can be specified with the flag
``-r /usr/local/share/igblast/fasta/imgt_human_tr_v.fasta``.

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
    docker run -v $DATA_DIR:/data\:z immcantation/suite:|docker-version| \\
        presto-abseq -1 $READS_R1 -2 $READS_R2 -y $YAML \\
        -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation_suite-|docker-version|.sif \\
        presto-abseq -1 $READS_R1 -2 $READS_R2 -y $YAML \\
        -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC


.. _ClontechPipeline:

Takara Bio / Clontech SMARTer v1 immune sequencing kit preprocessing pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A start to finish pRESTO processing script for Takara Bio / Clontech SMARTer v1 immune
sequencing kit data. C-regions are assigned using the universal C-region primer sequences are
available from the Immcantation repository under
`protocols/Universal <https://bitbucket.org/kleinstein/immcantation/src/master/protocols/Universal>`__
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
    docker run -v $DATA_DIR:/data\:z immcantation/suite:|docker-version| \\
        presto-clontech -1 $READS_R1 -2 $READS_R2 -j $CREGION -r $VREF \\
        -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation_suite-|docker-version|.sif \\
        presto-clontech -1 $READS_R1 -2 $READS_R2 -j $CREGION -r $VREF \\
        -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC


.. _ClontechPipelineUMI:

Takara Bio / Clontech SMARTer v2 (UMI) immune sequencing kit preprocessing pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A start to finish pRESTO processing script for Takara Bio / Clontech SMARTer v2 immune
sequencing kit data that includes UMIs. C-regions are assigned using the universal C-region
primer sequences are available from the Immcantation repository under
`protocols/Universal <https://bitbucket.org/kleinstein/immcantation/src/master/protocols/Universal>`__
or inside the container under ``/usr/local/share/protocols/Universal``.

.. include:: ../_include/usage.rst
    :start-after: Start presto-clontech-umi
    :end-before: End presto-clontech-umi

**Example: presto-clontech-umi**

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
    docker run -v $DATA_DIR:/data\:z immcantation/suite:|docker-version| \\
        presto-clontech-umi -1 $READS_R1 -2 $READS_R2 -j $CREGION -r $VREF \\
        -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation_suite-|docker-version|.sif \\
        presto-clontech-umi -1 $READS_R1 -2 $READS_R2 -j $CREGION -r $VREF \\
        -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

.. _10XPipeline:

10x Genomics V(D)J annotation pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assigns new annotations and infers clonal relationships to 10x Genomics
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
    docker run -v $DATA_DIR:/data\:z immcantation/suite:|docker-version| \\
        changeo-10x -s $READS -a $ANNOTATIONS -x $DIST -n $SAMPLE_NAME \\
        -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation_suite-|docker-version|.sif \\
        changeo-10x -s $READS -a $ANNOTATIONS -x $DIST -n $SAMPLE_NAME \\
        -o $OUT_DIR -p $NPROC


.. _IgBLASTPipeline:

IgBLAST annotation pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    docker run -v $DATA_DIR:/data\:z immcantation/suite:|docker-version| \\
        changeo-igblast -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation_suite-|docker-version|.sif \\
        changeo-igblast -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC


.. _GenotypePipeline:

Genotyping pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    docker run -v $DATA_DIR:/data\:z immcantation/suite:|docker-version| \\
        tigger-genotype -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation_suite-|docker-version|.sif \\
        tigger-genotype -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC


TIgGER infers the subject-specific genotyped V gene calls and saves the corrected calls in a new column, ``v_call_genotyped``. 
TIgGER also generates a ``*_genotype.fasta`` file, which contains the subject-specific germline IGHV genes. In future analyses, 
if ``v_call_genotyped`` column is used to replace ``v_call``, please remember to use  this ``*_genotype.fasta`` file generated previously 
by TIgGER as the subject-specific IGHV gene germline. An example of this application can be found in the **Clonal assignment pipeline** section. 


.. _ThresholdPipeline:

Clonal threshold inference pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    docker run -v $DATA_DIR:/data\:z immcantation/suite:|docker-version| \\
        shazam-threshold -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation_suite-|docker-version|.sif \\
        shazam-threshold -d $DB -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC


.. _ClonePipeline:

Clonal assignment pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assigns Ig sequences into clonally related lineages and builds full germline
sequences.

If the TIgGER, or another package, was applied previously to the data set for 
identifying a subject-specific genotype, including potentially novel V, D 
and/or J genes, a new directory $NEW_REF with the personalized germline database 
should be created.  For example, if TIgGER was run to identify a subject-specific 
IGHV genotype, the directory would contain: 1) ``*_genotype.fasta`` file generated 
previously by TIgGER, which contains the subject-specific germline IGHV genes 
2) ``imgt_human_IGHD.fasta`` and ``imgt_human _IGHJ.fasta``, which contain the IMGT IGHD 
and IGHJ genes and can both be copied from the original germline 
database: ``/usr/local/share/germlines/imgt/human/vdj/``. When changeo-clone is called, 
this new personalized germline database should be passed with parameter ``-r`` 
(see example below). And please remember to update ``v_call`` column with 
subject-specific IGHV call (for TIgGER this is found in ``v_call_genotyped`` column).

.. code-block:: R

   # update v_call
   db %>% 
      dplyr::mutate(v_call = v_call_genotyped) %>% 
      select(-v_call_genotyped)  

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
    docker run -v $DATA_DIR:/data\:z immcantation/suite:|docker-version| \\
        changeo-clone -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation_suite-|docker-version|.sif \\
        changeo-clone -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

**Example: changeo-clone with personalized germline database**

.. parsed-literal::

    # Arguments
    DATA_DIR=~/project
    NEW_REF=/data/personalized_germlines
    DB=/data/changeo/sample/sample_genotyped.tab
    DIST=0.15
    SAMPLE_NAME=sample
    OUT_DIR=/data/changeo/sample
    NPROC=4

    # Run pipeline in docker image
    docker run -v $DATA_DIR:/data:z immcantation/suite:|docker-version| \
    changeo-clone -r $NEW_REF -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

    # Singularity command
    singularity exec -B $DATA_DIR:/data immcantation_suite-|docker-version|.sif \
    changeo-clone -r $NEW_REF -d $DB -x $DIST -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC
