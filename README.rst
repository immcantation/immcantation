Immcantation
============

Advances in high-throughput sequencing technologies now allow for
large-scale characterization of B cell receptor (BCR) and T cell
receptor (TCR) repertoires. The high germline and somatic diversity of
the adaptive immune receptor repertoire (AIRR) presents challenges for
biologically meaningful analysis - requiring the development of
specialized computational methods.

The Immcantation framework provide a start-to-finish analytical
ecosystem for high-throughput AIRR-seq datasets. Beginning from raw
reads, Python and R packages are provided for pre-processing, population
structure determination, and repertoire analysis.

Repository
~~~~~~~~~~

This repository contains common documentation, accessory scripts,
template pipelines, and docker build files for tools in the Immcantation
framework.

+-------------+-------------------------------------------------------------------------------------------------------+
| Folder      | Contents                                                                                              |
+=============+=======================================================================================================+
| docker      | Dockerfiles for images hosted on `Docker Hub <https://hub.docker.com/r/kleinstein/immcantation>`__.   |
+-------------+-------------------------------------------------------------------------------------------------------+
| docs        | Sphinx build files for docs hosted on `ReadTheDocs <https://immcantation.readthedocs.io>`__.          |
+-------------+-------------------------------------------------------------------------------------------------------+
| pipelines   | Pipeline template scripts for the docker images.                                                      |
+-------------+-------------------------------------------------------------------------------------------------------+
| protocols   | Primer sequences and amplicon designs for published experimental protocols.                           |
+-------------+-------------------------------------------------------------------------------------------------------+
| scripts     | Accessory scripts for IMGT, IgBLAST and VDJTools.                                                     |
+-------------+-------------------------------------------------------------------------------------------------------+

Docker Container
~~~~~~~~~~~~~~~~

We have provided a complete installation of the Immcantation framework,
its dependencies, accessory scripts, and IgBLAST in a
`Docker <http://www.docker.com>`__ image. The image also includes both
the IgBLAST and IMGT reference germline sets, as well as several
template pipeline scripts. The image is available on docker hub at
`kleinstein/immcantation <https://hub.docker.com/r/kleinstein/immcantation>`__

Images are versioned through tags with images containing official
releases denoted by meta-version numbers (eg, ``1.0.0``). The ``devel``
tag denotes the latest development (unstabled) builds.

Documentation
~~~~~~~~~~~~~

Complete usage documentation, API documentation, and several tutorials
for the Immcantation framework tools can be found on the `Immcantation
Portal <https://immcantation.readthedocs.io>`__.
