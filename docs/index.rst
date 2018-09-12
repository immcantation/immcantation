.. Immcantation: An Integrated Framework for Adaptive Immune Receptor Repertoire Analysis
.. keywords:  rep-seq, immuno-seq, repertoire sequencing, BCR, TCR, Ig, IgOme, AIRR, adaptive immunity,
    somatic hypermutation, vdj-seq, AbSeq, AbPair, iRepertoire, VDJ, Immunoglobulin

.. toctree::
    :maxdepth: 1
    :hidden:

    Welcome <self>
    Getting Started <intro>
    Contact & Cite <about>
    Contributing <contrib>

.. toctree::
    :maxdepth: 3
    :hidden:
    :caption: Docker Container

    docker/intro
    docker/guide
    docker/news

.. toctree::
    :maxdepth: 1
    :caption: Packages
    :hidden:

    pRESTO <http://presto.readthedocs.io>
    Change-O <http://changeo.readthedocs.io>
    Alakazam <http://alakazam.readthedocs.io>
    SHazaM <http://shazam.readthedocs.io>
    TIgGER <http://tigger.readthedocs.io>
    RDI <http://rdi.readthedocs.io>

.. toctree::
    :maxdepth: 1
    :caption: In Development
    :hidden:

    SCOPe <packages/scope>
    prestoR <packages/prestor>

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


**Click on the images below for more details.**

.. |presto| image:: _static/presto.png
    :align: middle
    :width: 200
    :target: http://presto.readthedocs.io

.. |changeo| image:: _static/changeo.png
    :align: middle
    :width: 200
    :target: http://changeo.readthedocs.io

.. |alakazam| image:: _static/alakazam.png
    :align: middle
    :width: 200
    :target: http://alakazam.readthedocs.io

.. |shazam| image:: _static/shazam.png
    :align: middle
    :width: 200
    :target: http://shazam.readthedocs.io

.. |tigger| image:: _static/tigger.png
    :align: middle
    :width: 200
    :target: http://tigger.readthedocs.io

.. |rdi| image:: _static/rdi.png
    :align: middle
    :width: 200
    :target: http://rdi.readthedocs.io

.. |scope| image:: _static/scope.png
    :align: middle
    :width: 200
    :target: packages/scope

.. |prestoR| image:: _static/prestoR.png
    :align: middle
    :width: 200
    :target: packages/prestor

.. list-table::
   :widths: 50 50

   * - |presto|
     - + Quality control
       + Read assembly
       + UMI processing
   * - |changeo|
     - + V(D)J reference alignment standardization
       + Clonal clustering
       + Conversion and annotation
   * - |alakazam|
     - + Clonal lineage reconstruction
       + Repertoire diversity
       + V(D)J gene usage analysis
       + Physicochemical property analysis
   * - |shazam|
     - + Mutation profiling
       + Selection analysis
       + SHM modeling
       + Clonal clustering threshold tuning
   * - |tigger|
     - + Novel polymorphism detection
       + Genotyping
   * - |rdi|
     - + Repertoire Dissimilarity Index
   * - |scope|
     - + Clonal clustering
   * - |prestoR|
     - + pRESTO report generation
