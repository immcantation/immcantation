Getting Started
===========================================================================================

Download & Installation
-------------------------------------------------------------------------------------------

Installation for both the Python and R packages is performed in the usual manner.

To install pRESTO and Change-O from PyPI::

    > pip3 install presto changeo --user

To install Alakazam, SHazaM, TIgGER and SCOPer from CRAN::

    > R
    > install.packages(c("alakazam", "shazam", "tigger", "scoper"))
    
Alternatively, a complete installation of the Immcantation framework and its dependencies
is available as a Docker container. Installation of the container is described
in :ref:`DockerIntro` and basic usage is described in :ref:`DockerGuide`.

Immcantation Tutorials
-------------------------------------------------------------------------------------------

Each tool in the framework has its own documentation site, with detailed usage information 
and examples. A good starting point to familiarize yourself with the framework is to
follow one the tutorials listed here.


Introductory Webinar and Jupyter Notebook
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a detailed use example for each Immcantation tool see the
`slides and example data <https://goo.gl/FpW3Sc>`_ from our introductory webinar series. 
You can find a Jupyter notebook version of the webinar `in the repository <https://bitbucket.org/kleinstein/immcantation/src/default/training/>`_. If you don't want to execute the Jupyter notebook yourself, you can explore a `website version of it here <https://kleinstein.bitbucket.io/tutorials/intro-lab/index.html>`_.

This webinar covers:

* V(D)J gene annotation and novel polymorphism detection

* Inference of B cell clonal relationships

* Diversity analysis

* Mutational load profiling

* Modeling of somatic hypermutation (SHM) targeting

* Quantification of selection pressure

10X Genomics BCR Repertoire Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For information on how to process 10X data to be analyzed with Immcantation, visit our :ref:`10X Genomics V(D)J sequence analysis tutorial<10X-walkthrough>` or the `Parsing 10X Genomics V(D)J data Change-O example page <https://changeo.readthedocs.io/en/stable/examples/10x.html>`_.


Overview of B Cell Repertoire Analysis
-------------------------------------------------------------------------------------------

    **Yaari and Kleinstein.**
    Practical guidelines for B-cell receptor repertoire sequencing analysis.
    *Genome Medicine. 7, 121 (2015).*
    `doi\:10.1186/s13073-015-0243-2 <http://doi.org/10.1186/s13073-015-0243-2>`__


