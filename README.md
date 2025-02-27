[![Docker Pulls](https://img.shields.io/docker/pulls/immcantation/suite)](https://hub.docker.com/u/immcantation)
[![](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)

# Immcantation

Advances in high-throughput sequencing technologies now allow for
large-scale characterization of B cell receptor (BCR) and T cell
receptor (TCR) repertoires. The high germline and somatic diversity of
the adaptive immune receptor repertoire (AIRR) presents challenges
for biologically meaningful analysis - requiring the development of
specialized computational methods.

The Immcantation framework provide a start-to-finish analytical
ecosystem for high-throughput AIRR-seq datasets. Although Immcantation
is focused on BCRs, methods are applicable to TCRs. Beginning from raw
reads, Python and R packages are provided for pre-processing,
population structure determination, and repertoire analysis.

### Repository

**IMPORTANT!** 
Immcantation has moved to https://github.com/immcantation/immcantation

To update Git configuration settings use:

```
   git config user.email "your-gh-user@email.com"
   git config user.name "your-gh-user-name"
   git remote set-url origin git@github.com:immcantation/immcantation.git
```

This repository contains common documentation, accessory scripts,
example pipelines, and docker build files for tools in the Immcantation
framework.

Folder      | Contents
----------- | ----------------------------------------------------------
docker      | Dockerfiles for images hosted on [Docker Hub](https://hub.docker.com/r/immcantation).
docs        | Sphinx build files for docs hosted on [ReadTheDocs](https://immcantation.readthedocs.io).
pipelines   | Example pipeline scripts for the docker images.
protocols   | Primer sequences and amplicon designs for published experimental protocols.
scripts     | Accessory scripts for IMGT, IgBLAST and VDJTools.

### Docker Container

We have provided a complete installation of the Immcantation framework,
its dependencies, accessory scripts, and IgBLAST in a
[Docker](http://www.docker.com) image. The image also includes both the
IgBLAST and IMGT reference germline sets, as well as several example
pipeline scripts. The image is available on Docker Hub at
[immcantation/suite](https://hub.docker.com/r/immcantation/suite)

Images are versioned through tags with images containing official
releases denoted by meta-version numbers (eg, `1.0.0`). The `devel` tag
denotes the latest development (unstable) builds. The tag `latest` is not
available. Images can be obtained with the command
`docker pull immcantation/suite:<tag>`.

### Documentation

Complete usage documentation, API documentation, and several tutorials
for the Immcantation framework tools can be found on the
[Immcantation Portal](https://immcantation.readthedocs.io).
