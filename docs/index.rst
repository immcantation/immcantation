.. meta::
   :description: Immcantation: An Integrated Framework for Adaptive Immune Receptor Repertoire Analysis
   :keywords: rep-seq, immuno-seq, vdj-seq, repertoire sequencing, BCR, TCR, Ig, AIRR,
    adaptive immunity, somatic hypermutation, AbSeq, AbPair, VDJ, immunoglobulin

.. meta::
    :twitter\:card:  summary_large_image
    :og\:title: Immcantation: An Integrated Framework for Adaptive Immune Receptor Repertoire Analysis
    :og\:image: _static/immcantation-card.png'

.. toctree::
    :maxdepth: 1
    :hidden:

    Welcome <self>
    Getting Started <intro>
    Tutorials <tutorials/tutorials>
    Data Standards <datastandards>
    Contact & Cite <about>
    Contributing <contrib>

.. toctree::
    :maxdepth: 3
    :hidden:
    :caption: Docker Container

    docker/intro
    docker/guide
    docker/pipelines
    docker/news

.. toctree::
    :maxdepth: 1
    :caption: Core Packages
    :hidden:

    pRESTO <http://presto.readthedocs.io>
    Change-O <http://changeo.readthedocs.io>
    Alakazam <http://alakazam.readthedocs.io>
    SHazaM <http://shazam.readthedocs.io>
    TIgGER <http://tigger.readthedocs.io>
    SCOPer <http://scoper.readthedocs.io>
    dowser <https://dowser.readthedocs.io/>

.. toctree::
    :maxdepth: 1
    :caption: Contributed Packages
    :hidden:

    RDI <http://rdi.readthedocs.io>
    RAbHIT <https://yaarilab.bitbucket.io/RAbHIT/>
    IgPhyML <https://igphyml.readthedocs.io>
    sumrep <https://github.com/matsengrp/sumrep>

.. toctree::
    :maxdepth: 1
    :caption: In Development
    :hidden:

    prestoR <packages/prestor>
    enchantR <https://enchantr.readthedocs.io/>


.. _Welcome:

Welcome to the Immcantation Portal!
==========================================================================================

Advances in high-throughput sequencing technologies now allow for large-scale
characterization of B cell receptor (BCR) and T cell receptor (TCR) repertoires. The high
germline and somatic diversity of the adaptive immune receptor repertoire (AIRR) presents
challenges for biologically meaningful analysis - requiring the development of specialized
computational methods.

The Immcantation framework provide a start-to-finish analytical ecosystem for
high-throughput AIRR-seq datasets. Beginning from raw reads, Python and R packages are
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

The different tools are available from PyPi, CRAN and Bitbucket. Versioned containers with
all tools installed are hosted on `Docker Hub <https://hub.docker.com/r/immcantation/suite>`_.


Core Packages
-----------------------------------------------------------------------------------------

**Click on the images below for more details.**

.. list-table::
    :widths: 40 60
    :align: left

    * - |presto-img|
      - .. image:: https://img.shields.io/pypi/dm/presto
            :target: https://pypi.org/project/presto
            :align: right
            :alt: downloads

        **pRESTO**

        + Quality control
        + Read assembly
        + UMI processing
        + Error profiling

    * - |changeo-img|
      - .. image:: https://img.shields.io/pypi/dm/changeo
            :target: https://pypi.org/project/changeo
            :align: right
            :alt: downloads

        **Change-O**

        + V(D)J reference alignment standardization
        + Clonal clustering
        + Germline reconstruction
        + Conversion and annotation

    * - |alakazam-img|
      - .. image:: https://cranlogs.r-pkg.org/badges/alakazam
            :target: https://www.r-pkg.org/pkg/alakazam
            :align: right
            :alt: downloads

        **Alakazam**

        + Clonal lineage reconstruction
        + Lineage topology analysis
        + Repertoire diversity
        + V(D)J gene usage
        + Physicochemical property analysis

    * - |shazam-img|
      - .. image:: https://cranlogs.r-pkg.org/badges/shazam
            :target: https://www.r-pkg.org/pkg/shazam
            :align: right
            :alt: downloads

        **SHazaM**

        + Mutation profiling
        + Selection pressure quantification
        + Empirical SHM models
        + Chimera detection
        + Clonal clustering threshold tuning

    * - |tigger-img|
      - .. image:: https://cranlogs.r-pkg.org/badges/tigger
            :target: https://www.r-pkg.org/pkg/tigger
            :align: right
            :alt: downloads

        **TIgGER**

        + Novel polymorphism detection
        + Genotyping

    * - |scoper-img|
      - .. image:: https://cranlogs.r-pkg.org/badges/scoper
            :target: https://www.r-pkg.org/pkg/scoper
            :align: right
            :alt: downloads

        **SCOPer**

        + Spectral clonal clustering methods

    * - |prestoR-img|
      - **prestoR**

        + pRESTO report generation

    * - |dowser-img|
      - .. image:: https://cranlogs.r-pkg.org/badges/dowser
            :target: https://www.r-pkg.org/pkg/dowser
            :align: right
            :alt: downloads

        **dowser**

        + Phylogenetic analysis on B cell receptor repertoires


Contributed Packages
-----------------------------------------------------------------------------------------

**Click on the images below for more details.**

.. list-table::
    :widths: 40 60
    :align: left

    * - |rdi-img|
      - **RDI**

        + Repertoire Dissimilarity Index

    * - |rabhit-img|
      - **RAbHIT**

        + Determination of V-D-J haplotypes

    * - |igphyml-img|
      - **IgPhyML**

        + Clonal lineage tree construction
        + Mutation/selection hypothesis testing

    * - |sumrep-img|
      - **sumrep**

        + Generate repertoire summary statistics.
        + Visualize and comparing repertoire summaries.

.. Image substitutions

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

.. |sumrep-img| image:: _static/sumrep.png
    :align: middle
    :width: 180
    :target: sumrep_
    :alt: sumrep

.. Doc links

.. _pRESTO: https://presto.readthedocs.io
.. _Change-O: https://changeo.readthedocs.io
.. _Alakazam: https://alakazam.readthedocs.io
.. _SHazaM: https://shazam.readthedocs.io
.. _TIgGER: https://tigger.readthedocs.io
.. _SCOPer: https://scoper.readthedocs.io
.. _dowser: https://dowser.readthedocs.io
.. _RDI: https://rdi.readthedocs.io
.. _IgPhyML: https://igphyml.readthedocs.io
.. _RAbHIT: https://yaarilab.bitbucket.io/RAbHIT/
.. _sumrep: https://github.com/matsengrp/sumrep
.. _prestoR: packages/prestor.html
