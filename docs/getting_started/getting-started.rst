Immcantation Tutorials
===========================================================================================

Each tool in the framework has its own documentation site, with detailed usage information
and examples. A good starting point to familiarize yourself with the framework is to
follow one the tutorials listed here. They are available as Notebooks under the
folder `training <https://bitbucket.org/kleinstein/immcantation/src/master/training>`_
in the Immcantation Bitbucket repository.

The Immcantation Lab container
-------------------------------------------------------------------------------------------

.. image:: ../_static/immcantationlab.png
    :alt: How to use the Immcantation Lab container
    :align: left

Information on obtaining and using the container with the notebooks and data needed to follow 
the Immcantation tutorials can be found within the general :ref:`guide <DockerGuide>`.


Introductory Webinar
-------------------------------------------------------------------------------------------

.. image:: ../_static/start-here.png
    :alt: Start here: Introductory Webinar and R Notebook
    :align: left

For a detailed use example for each Immcantation tool follow the tutorial
*Introduction to B cell repertoire analysis using the Immcantation
framework*. This tutorial is based on our
introductory webinar. It is also available as an R Notebook
(`intro-lab.Rmd <https://bitbucket.org/kleinstein/immcantation/src/master/training>`_) in the Immcantation repository.

.. admonition:: Tutorial

    :doc:`Introduction to B cell repertoire analysis using the Immcantation framework <intro-lab>` covers:

    + V(D)J gene annotation and novel polymorphism detection
    + Inference of B cell clonal relationships
    + Diversity analysis
    + Mutational load profiling
    + Modeling of somatic hypermutation (SHM) targeting
    + Quantification of selection pressure



Single-cell Analysis
-------------------------------------------------------------------------------------------


.. image:: ../_static/bcell.png
    :alt: 10x Genomics V(D)J Sequence Analysis with Immcantation Tutorial
    :align: left

For information on how to process 10x Genomics VDJ data to be analyzed with Immcantation, we
offer the introductory tutorial *10x Genomics V(D)J Sequence Analysis with Immcantation*. It is 
available as an R Notebook (`10x_tutorial.Rmd <https://bitbucket.org/kleinstein/immcantation/src/master/training>`_) 
in the Immcantation repository.

.. admonition:: Tutorial

    :doc:`10x Genomics V(D)J Sequence Analysis with Immcantation <10x_tutorial>` covers:

    + V(D)J gene annotation
    + Inference of clonal relationships
    + V gene Somatic Hypermutation (SHM) frequency
    + Lineage tree reconstruction
    + Incorporation of Cell Ranger annotations

.. image:: ../_static/bcellgex.png
    :alt: Integration of BCR and GEX data
    :align: left

In *Integration of BCR and GEX data* we demonstrate an enhanced analysis by integrating 10x BCR and 10x GEX data. 
The R Notebook (`BCR_Seurat_tutorial.Rmd <https://bitbucket.org/kleinstein/immcantation/src/master/training/>`_) 
is available in the Immcantation repository.

.. admonition:: Tutorial

    :doc:`Integration of BCR and GEX data <BCR_Seurat_tutorial>` covers:

    + Integration of BCR data with the GEX Seurat object
    + Highlight BCR cells in the GEX UMAP
    + Integration of GEX cell annotations in the BCR data
    + Identify GEX clusters in the BCR UMAP
    + Highlight other BCR features in UMAPs


Video presentations
===========================================================================================

.. image:: ../_static/immcantation-yt.png
    :target: https://www.youtube.com/channel/UCWQgFmSnv8B9Q5G_kunOabw/playlists
    :alt: Link to Immcantation's YouTube Channel
    :align: left

You can watch presentations by Immcantation developers and users in Immcantation's YouTube channel.


Vignettes
===========================================================================================

Detailed usage documentation and tutorials for each individual tool in Immcantation are
provided in the main documentation pages for each tool. The following list of shortcuts
cover common analyses. Note, each link will leave the Immcantation portal page.

.. toctree::
    :maxdepth: 1

    Assembling raw reads from simple Illumina sequencing protocols with pRESTO <https://presto.readthedocs.io/en/stable/workflows/Greiff2014_Workflow.html>
    Assembling raw reads from 5'RACE UMI barcoded Illumina sequencing protocols with pRESTO <https://presto.readthedocs.io/en/stable/workflows/VanderHeiden2017_Workflow.html>
    Processing 10x Genomics Cell Ranger data with Change-O <https://changeo.readthedocs.io/en/stable/examples/10x.html>
    Processing IgBLAST data with Change-O <https://changeo.readthedocs.io/en/stable/examples/igblast.html>
    Processing IMGT/HighV-QUEST data with Change-O <https://changeo.readthedocs.io/en/stable/examples/imgt.html>
    Building lineage trees with IgPhyML <https://changeo.readthedocs.io/en/stable/examples/igphyml.html>
    Assigning clonal groups with SCOPer <https://scoper.readthedocs.io/en/stable/vignettes/Scoper-Vignette>
    Basic gene usage analysis with Alakazam <https://alakazam.readthedocs.io/en/stable/vignettes/GeneUsage-Vignette>
    Clonality and diversity analysis with Alakazam <https://alakazam.readthedocs.io/en/stable/vignettes/Diversity-Vignette>
    Mutational load analysis with SHazaM <https://shazam.readthedocs.io/en/stable/vignettes/Mutation-Vignette>
    Selection pressure analysis with SHazaM <https://shazam.readthedocs.io/en/stable/vignettes/Baseline-Vignette>
    Building SHM targeting models with SHazaM <https://shazam.readthedocs.io/en/stable/vignettes/Targeting-Vignette>
    Novel allele detection and genotyping with TIgGER <https://tigger.readthedocs.io/en/stable/vignettes/Tigger-Vignette>
