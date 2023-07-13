.. _DockerIntro:

Docker Container Installation and Overview
================================================================================

We have provided a complete installation of the Immcantation framework, its
dependencies, accessory scripts, and IgBLAST in a
`Docker container <http://www.docker.com>`__. The `container <https://hub.docker.com/r/immcantation/suite/>`__ also includes both the IgBLAST and
IMGT reference germline sets, as well as several example pipeline scripts.

We currently have four containers available on `DockerHub <https://hub.docker.com/r/immcantation/>`__:

+---------------------------------------+-----------------------------------------------------------------------------------------+
| Name                                  | Contents                                                                                |
+=======================================+=========================================================================================+
| immcantation/suite                    | Immcantation suite, supporting applications and databases.                              |
+---------------------------------------+-----------------------------------------------------------------------------------------+
| immcantation/lab                      | Immcantation tutorial materials. Only for training, not to be used in production.       |
+---------------------------------------+-----------------------------------------------------------------------------------------+
| immcantation/base                     | Base image for Immcantation builds.                                                     |
+---------------------------------------+-----------------------------------------------------------------------------------------+
| immcantation/test                     | Immcantation unit test image.                                                           |
+---------------------------------------+-----------------------------------------------------------------------------------------+

**For tutorial purposes**, use immcantation/lab (be sure to replace ``suite`` with ``lab`` in the following code chunks) and follow the directions :ref:`here <DockerGuideTutorials>`. For all Immcantation uses, use immcantation/suite.
Note that containers are versioned through tags with containers containing official releases
denoted by meta-version numbers (``x.y.z``). The ``devel`` tag denotes the
latest development (unstable) builds.

.. include:: installation.rst

.. include:: overview.rst

.. include:: guide.rst

.. include:: pipelines.rst

.. include:: news.rst
