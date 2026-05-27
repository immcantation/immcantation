.. meta::
   :description: Immcantation: An Integrated Framework for Adaptive Immune Receptor Repertoire Analysis
   :keywords: AIRR-seq, rep-seq, immuno-seq, vdj-seq, repertoire sequencing, BCR, TCR, Ig, AIRR,
    adaptive immunity, somatic hypermutation, SHM, AbSeq, AbPair, VDJ, immunoglobulin, bulk, single cell, scRNA-seq

.. meta::
    :twitter\:card:  summary_large_image
    :og\:title: Immcantation: An Integrated Framework for Adaptive Immune Receptor Repertoire Analysis
    :og\:image: _static/immcantation-card.png

.. toctree::
    :maxdepth: 1
    :hidden:

    Welcome <self>
    Data Standards <datastandards>
    Contact & Cite <about>
    Contributing <contrib>

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Installation

    intro
    docker/intro
    docker/news

.. toctree::
    :maxdepth: 3
    :hidden:
    :caption: Getting started

    getting_started/getting-started
    getting_started/resources

.. toctree::
    :maxdepth: 1
    :caption: Core Packages
    :hidden:

    Alakazam <https://alakazam.readthedocs.io>
    Change-O <https://changeo.readthedocs.io>
    Dowser <https://dowser.readthedocs.io/>
    pRESTO <https://presto.readthedocs.io>
    SCOPer <https://scoper.readthedocs.io>
    SHazaM <https://shazam.readthedocs.io>
    TIgGER <https://tigger.readthedocs.io>

.. toctree::
    :maxdepth: 1
    :caption: Contributed Packages
    :hidden:

    IgPhyML <https://igphyml.readthedocs.io>
    PIgLET <https://bitbucket.org/yaarilab/piglet>
    RAbHIT <https://yaarilab.bitbucket.io/RAbHIT/>
    RDI <https://rdi.readthedocs.io>
    sumrep <https://github.com/matsengrp/sumrep>

.. toctree::
    :maxdepth: 1
    :caption: Workflows
    :hidden:

    nf-core/airrflow <https://nf-co.re/airrflow>

.. toctree::
    :maxdepth: 1
    :caption: In Development
    :hidden:

    enchantR <https://enchantr.readthedocs.io/>
    prestoR <packages/prestor>


.. _Welcome:

Welcome to the Immcantation Portal!
==========================================================================================

Advances in high-throughput sequencing technologies now allow for large-scale
characterization of B cell receptor (BCR) and T cell receptor (TCR) repertoires. The high
germline and somatic diversity of the adaptive immune receptor repertoire (AIRR) presents
challenges for biologically meaningful analysis - requiring the development of specialized
computational methods.

The Immcantation framework provides a **start-to-finish analytical ecosystem for
high-throughput AIRR-seq** datasets. Beginning from raw reads, Python and R packages are
provided for pre-processing, population structure determination, and repertoire analysis.

.. image:: https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic
    :target: https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html
    :align: left

Immcantation supports both the original Change-O standard and the new Adaptive Immune
Receptor Repertoire (AIRR) standard developed by the
`AIRR Community (AIRR-C) <https://www.antibodysociety.org/the-airr-community/>`_.

.. image:: https://img.shields.io/docker/pulls/immcantation/suite
    :target: https://hub.docker.com/u/immcantation
    :align: left

The different tools are available from PyPi, CRAN and GitHub. Versioned containers with
all tools installed are hosted on `Docker Hub <https://hub.docker.com/r/immcantation/suite>`_.
A best practices workflow using the Immcantation packages is available as 
a `nf-core/airrflow <https://nf-co.re/airrflow>`_ Nextflow workflow.


Start-to-finish Workflow
------------------------------------------------------------------------------------------
Run the Immcantation packages in a single workflow to analyze your
repertoire sequencing data.

**Click on the image below for more details.**

.. list-table::
    :widths: 40 60
    :align: left



    * - |airrflow-img|
      
      - **nf-core/airrflow**

        + Nextflow workflow using Immcantation
        + Bulk and single-cell BCR and TCR analysis
        + Part of the `nf-core <https://nf-co.re/>`_ project

.. raw:: html

   <div style="clear: left;"></div>

Core Packages
-----------------------------------------------------------------------------------------

Core packages are the software tools maintained by the Immcantation team. 
These packages provide essential functionality for processing, analyzing, and visualizing 
AIRR-seq data.

**Click on the images below for more details.**

.. list-table::
    :widths: 40 60
    :align: left



    * - |alakazam-img|
      - .. image:: https://cranlogs.r-pkg.org/badges/alakazam
            :target: https://www.r-pkg.org/pkg/alakazam
            :align: right
            :alt: downloads

        **Alakazam**

        + Repertoire diversity
        + V(D)J gene usage
        + Physicochemical property analysis

        .. raw:: html

           <div style="text-align: right;"><a href="https://medicine.yale.edu/lab/kleinstein/">Kleinstein Lab</a></div>

    * - |amulety-img|
      - .. image:: https://img.shields.io/pypi/dm/amulety
            :target: https://pypi.org/project/amulety
            :align: right
            :alt: downloads

        **Amulety**

        + BCR and TCR embeddings
        + Multiple embedding models supported

        .. raw:: html

           <div style="text-align: right;"><a href="https://medicine.yale.edu/lab/kleinstein/">Kleinstein Lab</a></div>

    * - |changeo-img|
      - .. image:: https://img.shields.io/pypi/dm/changeo
            :target: https://pypi.org/project/changeo
            :align: right
            :alt: downloads

        **Change-O**

        + V(D)J alignment with IgBLAST and IMGT

        .. raw:: html

           <div style="text-align: right;"><a href="https://medicine.yale.edu/lab/kleinstein/">Kleinstein Lab</a></div>

    * - |dowser-img|
      - .. image:: https://cranlogs.r-pkg.org/badges/dowser
            :target: https://www.r-pkg.org/pkg/dowser
            :align: right
            :alt: downloads

        **Dowser**

        + B cell lineage trees
        + Migration and differentiation analysis
        + Detect ongoing evolution over time

        .. raw:: html

           <div style="text-align: right;"><a href="https://sites.dartmouth.edu/hoehn/">Hoehn lab</a></div>

    * - |presto-img|
      - .. image:: https://img.shields.io/pypi/dm/presto
            :target: https://pypi.org/project/presto
            :align: right
            :alt: downloads

        **pRESTO**

        + Bulk BCR sequence data pre-processing
        + Read assembly and QC
        + UMI processing

        .. raw:: html

           <div style="text-align: right;"><a href="https://medicine.yale.edu/lab/kleinstein/">Kleinstein Lab</a></div>

    * - |scoper-img|
      - .. image:: https://cranlogs.r-pkg.org/badges/scoper
            :target: https://www.r-pkg.org/pkg/scoper
            :align: right
            :alt: downloads

        **SCOPer**

        + Identify clonal relationships

        .. raw:: html

           <div style="text-align: right;"><a href="https://medicine.yale.edu/lab/kleinstein/">Kleinstein Lab</a></div>

    * - |shazam-img|
      - .. image:: https://cranlogs.r-pkg.org/badges/shazam
            :target: https://www.r-pkg.org/pkg/shazam
            :align: right
            :alt: downloads

        **SHazaM**

        + Clonal clustering threshold tuning
        + Mutation profiling
        + Selection pressure quantification
        + Empirical SHM models

        .. raw:: html

           <div style="text-align: right;"><a href="https://medicine.yale.edu/lab/kleinstein/">Kleinstein Lab</a></div>

    * - |tigger-img|
      - .. image:: https://cranlogs.r-pkg.org/badges/tigger
            :target: https://www.r-pkg.org/pkg/tigger
            :align: right
            :alt: downloads

        **TIgGER**

        + Novel polymorphism detection
        + Genotyping

        .. raw:: html

           <div style="text-align: right;"><a href="https://medicine.yale.edu/lab/kleinstein/">Kleinstein Lab</a></div>

Contributed Packages
-----------------------------------------------------------------------------------------

Contributed packages in the Immcantation ecosystem are immunoinformatics software packages 
developed, shared, and maintained by the community. These packages interoperate with the 
Immcantation framework through the `AIRR Community Standard <https://docs.airr-community.org/en/stable/datarep/rearrangements.html>`__ 
and complement Immcantation by providing specialized functionality for AIRR analysis.

**Click on the images below for more details.**

.. list-table::
    :widths: 40 60
    :align: left

    * - |igphyml-img|
      - **IgPhyML**

        + Method to build lineage trees
        + Mutation/selection hypothesis testing
        + Best used via `Dowser`_ package

        .. raw:: html

           <div style="text-align: right;"><a href="https://sites.dartmouth.edu/hoehn/">Hoehn Lab</a></div>

    * - |rabhit-img|
      - **RAbHIT**

        + Determination of V-D-J haplotypes

        .. raw:: html

           <div style="text-align: right;"><a href="https://medicine.yale.edu/profile/gur-yaari/">Yaari Lab</a></div>

    * - |piglet-img|
      - **PIgLET**

        + Tools to improve genotype inference      

        .. raw:: html

           <div style="text-align: right;"><a href="https://medicine.yale.edu/profile/gur-yaari/">Yaari Lab</a></div>

    * - |rdi-img|
      - **RDI**

        + Repertoire Dissimilarity Index

        .. raw:: html

           <div style="text-align: right;"><a href="https://medicine.yale.edu/lab/kleinstein/">Kleinstein Lab</a></div>

    * - |sumrep-img|
      - **sumrep**

        + Generate repertoire summary statistics.
        + Visualize and comparing repertoire summaries.

        .. raw:: html

           <div style="text-align: right;"><a href="https://matsen.fhcrc.org/">Matsen Group</a></div>

.. Image substitutions

.. |airrflow-img| image:: _static/airrflow_logo.png
    :align: middle
    :width: 200
    :target: airrflow_
    :alt: nf-core/airrflow

.. |presto-img| image:: _static/presto.png
    :align: middle
    :width: 200
    :target: pRESTO_
    :alt: pRESTO

.. |changeo-img| image:: _static/changeo.png
    :align: middle
    :width: 200
    :target: Change-O_
    :alt: Change-O

.. |alakazam-img| image:: _static/alakazam.png
    :align: middle
    :width: 200
    :target: Alakazam_
    :alt: alakazam

.. |amulety-img| image:: _static/amulety.png
    :align: middle
    :width: 200
    :target: Amulety_
    :alt: amulety

.. |shazam-img| image:: _static/shazam.png
    :align: middle
    :width: 200
    :target: SHazaM_
    :alt: SHazaM

.. |tigger-img| image:: _static/tigger.png
    :align: middle
    :width: 200
    :target: TIgGER_
    :alt: TIgGER

.. |scoper-img| image:: _static/scoper.png
    :align: middle
    :width: 200
    :target: SCOPer_
    :alt: SCOPer

.. |prestoR-img| image:: _static/prestoR.png
    :align: middle
    :width: 200
    :target: prestoR_
    :alt: prestoR

.. |dowser-img| image:: _static/dowser.png
    :align: middle
    :width: 200
    :target: dowser_
    :alt: dowser

.. |rdi-img| image:: _static/rdi.png
    :align: middle
    :width: 200
    :target: RDI_
    :alt: RDI

.. |igphyml-img| image:: _static/igphyml.png
    :align: middle
    :width: 180
    :target: IgPhyML_
    :alt: IgPhyML

.. |rabhit-img| image:: _static/rabhit.png
    :align: middle
    :width: 140
    :target: RAbHIT_
    :alt: RAbHIT

.. |piglet-img| image:: _static/piglet.png
    :align: middle
    :width: 140
    :target: PIgLET_
    :alt: PIgLET    

.. |sumrep-img| image:: _static/sumrep.png
    :align: middle
    :width: 180
    :target: sumrep_
    :alt: sumrep

.. Doc links

.. _airrflow: https://nf-co.re/airrflow
.. _Alakazam: https://alakazam.readthedocs.io
.. _Amulety: https://amulety.readthedocs.io
.. _Change-O: https://changeo.readthedocs.io
.. _Dowser: https://dowser.readthedocs.io
.. _IgPhyML: https://igphyml.readthedocs.io
.. _PIgLET: https://bitbucket.org/yaarilab/piglet
.. _pRESTO: https://presto.readthedocs.io
.. _prestoR: packages/prestor.html
.. _RAbHIT: https://yaarilab.bitbucket.io/RAbHIT/
.. _RDI: https://rdi.readthedocs.io
.. _SCOPer: https://scoper.readthedocs.io
.. _SHazaM: https://shazam.readthedocs.io
.. _sumrep: https://github.com/matsengrp/sumrep
.. _TIgGER: https://tigger.readthedocs.io
