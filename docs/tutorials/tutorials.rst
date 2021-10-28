Immcantation Tutorials
===========================================================================================

Each tool in the framework has its own documentation site, with detailed usage information 
and examples. A good starting point to familiarize yourself with the framework is to
follow one the tutorials listed here.


Introductory Webinar and Jupyter Notebook
-------------------------------------------------------------------------------------------

For a detailed use example for each Immcantation tool see the Jupyter notebook from our
introductory webinar `in the repository <https://bitbucket.org/kleinstein/immcantation/src/master/training/>`_.
If you don't want to execute the Jupyter notebook yourself, you can explore a
`website version of it here <https://kleinstein.bitbucket.io/tutorials/intro-lab/index.html>`_.
*Introduction to B cell repertoire analysis using the Immcantation framework* covers:

+ V(D)J gene annotation and novel polymorphism detection
+ Inference of B cell clonal relationships
+ Diversity analysis
+ Mutational load profiling
+ Modeling of somatic hypermutation (SHM) targeting
+ Quantification of selection pressure


Single-cell Analysis
-------------------------------------------------------------------------------------------

For information on how to process 10x Genomics VDJ data to be analyzed with Immcantation, we
offer an introductory tutorial for new users:

.. toctree::
    :maxdepth: 1

    10x_tutorial

We also demonstrate an enhanced analysis by integrating 10x BCR and 10x GEX data in 
this `jupyter notebook <https://bitbucket.org/kleinstein/immcantation/src/master/training/BCR_Seurat_tutorial.ipynb?viewer=nbviewer>`_.

Video presentations
-------------------------------------------------------------------------------------------

.. image:: ../_static/immcantation-yt.png
    :target: https://www.youtube.com/channel/UCWQgFmSnv8B9Q5G_kunOabw/about
    :alt: Link to Immcantation's YouTube Channel
    :align: left

You can watch presentations by Immcantation developers and users in Immcantation's YouTube channel.

Here we provide links to other publicly available Immcantation presentations:

`A unified maximum likelihood framework for B cell repertoire phylogenetics <https://www.youtube.com/watch?v=oIFDLP5UqLo&ab_channel=ISCB>`_
Kenneth B. Hoehn
ISMB 2018 (July 6-10, 2018)

`Repertoire-wide phylogenetic models of B-cell molecular evolution reveal evolutionary signatures of aging and vaccination <https://www.youtube.com/watch?v=GwXfXD84K7o&ab_channel=AIRRCommunity>`_
Kenneth B. Hoehn
AIRR Community Meeting IV (May 11-15, 2019)

`​B-cell lineage analysis and COVID-19 <https://www.youtube.com/watch?v=d1IRY0LSuv4&ab_channel=AIRRCommunity>`_
Kenneth B. Hoehn
AIRR Community Special Event: Leveraging AIRR-sequencing data to inform the biology of COVID-19 (September 8-10, 2020)

`Comprehensive analysis of immunogenomics sequencing data in the cloud to facilitate reproducibility and rigor of immunogenomics research <https://www.youtube.com/watch?v=lmR6NX4qXaw&ab_channel=ISCB>`_
Victor Greiff, Kenneth B. Hoehn, Steven H. Kleinstein, Serghei Mangul, Milena Pavlovic, Kerui Peng
Tutorials - ISMB/ECCB 2021 (July 22, 2021)


Vignettes
-------------------------------------------------------------------------------------------

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

