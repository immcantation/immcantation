# Immcantation training materials

[TOC]

Target audience: Bioinformaticians, immunologists, biologist, scientists and students learning about immune repertoires.

Prerequisite Skills and Knowledge Required: Familiarity with R, python and Linux.

Licence: These tutorials are licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License](https://creativecommons.org/licenses/by-nc/4.0/).

## Introduction to B cell repertoire analysis

Get a global overview of how the different tools in the Immcantation framework work together with a [Jupyter notebook](intro-lab.ipynb?viewer=nbviewer) based on the materials presented in the [webinar](https://immcantation.eventbrite.com). Use it online with
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/immcantation/immcantation-lab/master) or locally, following the instructions below (See [Using the modules locally with Docker](#markdown-header-using-the-modules-locally-with-docker)). Learn more about Jupyter notebooks [here](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/). If you don’t want to execute the Jupyter notebook at this time, you can explore [the website version, with the results, here](https://kleinstein.bitbucket.io/tutorials/intro-lab/index.html).

### Learning Outcomes

* V(D)J gene annotation and novel polymorphism detection
* Clonotype assignment
* Diversity analysis
* Mutational load profiling
* Modeling of somatic hypermutation targeting
* Quantification of selection pressure

### Domain Problem

The field of high-throughput adaptive immune receptor repertoire sequencing (AIRR-seq) has experienced significant growth in recent years, but this growth has come with considerable complexity and variety in experimental design. These complexities, combined with the high germline and somatic diversity of immunoglobulin repertoires, present analytical challenges requiring specialized methodologies. This tutorial will cover common investigative approaches and pitfalls in AIRR-seq data analysis.

### Dataset for the case study

Processed reads (input.fasta) from one healthy donor (PGP1) 3 weeks after flu vaccination (Laserson et al. (2014)). Download the introductory tutorial data [here](https://yale.box.com/shared/static/4bo611b70x8u92qvss1pypmcr9wmqil4).

### Funding

Developed with the support of the National Library of Medicine (NIH NLM T15 LM007056) and the National Institute of Allergy and Infectious Diseases (NIH NIAID R01 AI104739).

## 10x Genomics V(D)J analysis

Learn how to process 10x Genomics VDJ data to be analyzed with Immcantation.

Download the 10x Genomics V(D)J analysis data [here](http://clip.med.yale.edu/immcantation/examples/10x_data_2subj.zip).

## Integration of BCR data and GEX data

Learn how to integrate BCR repertoire and gene expression analysis with this [Jupyter notebook](BCR_Seurat_tutorial.ipynb?viewer=nbviewer).


### Learning Outcomes

* Integrate BCR data to GEX data in Seurat object
* Integrate annotations of cells from GEX data to BCR data

## Lineage tree reconstruction

Beginning with processed single cell RNA-seq (scRNA-seq) + BCR data from 10X Genomics, you will learn:

* how cell type annotations can be associated with BCR sequences,
* how clonal clusters can be identified, and
* how B cell phylogenetic trees can be built and visualized using these data sources.

# Using the tutorials

## With the Immcantation training Docker container

The training container `immcantation/lab`
contains Immcantation, Jupyter Notebook,
the training notebooks and all
required example data and software dependencies
ready to be used. This container is *not* to be
used in real analysis: the reference germline
IGHV3-20*04 has been removed to demonstrate how
TIgGER can identify this allele.
.

* 1 Pull the Immcantation Lab container image:

```
# Example: pull the development version of the lab
docker pull immcantation/lab:devel
```

* 2 Run the container:

```
docker run --network=host -it --rm -p 8888:8888 immcantation/lab:devel
```

Or, if you want to save the results in your computer:

```
# Note: change my-out-dir for the full path to the local directory where
# you want to have the results saved to
docker run --network=host -it --rm -v my-out-dir:/home/magus/notebooks/results:z -p 8888:8888 immcantation/lab:devel
```

Once the container is running, You will see in the terminal a message asking you to visit a url like `http://<hostname>:8888/?token=<token>`

* 3 Open your computer's internet browser and visit the url

When you visit the url from the previous step, you will start a Jupyter session in your browser.

```
# Example: http://localhost:8888/?token=18303237b2521e72f00685e4fdf754f955f82a958a8e57ec
```

* 4 Open the notebook

Open the notebook you want to work with. Use CTRL+Enter to execute the commands inside the cells.

For an introduction to Jupyter, visit the [official documentation site](https://jupyter-notebook.readthedocs.io/en/latest/).

## With a release container

The Immcantation release containers don't have Jupyter Notebook,
example data and tools used in the notebooks (like `Seurat`).
But you can download the example data you need, start the
release container, and copy/paste the commands from the tutorial into
the container terminal or R session.

```
## Get the 4.3.0 release of immcantation
# Note: For some operating systems, it may be necessary to use super-user
# privileges (sudo), and/or to have Docker Desktop running
docker pull immcantation/suite:4.3.0

## Run the container
# Note: this assumes you have the example data in the current directory

# Docker on Linux/Mac OS X
docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.3.0 bash

# Docker on Windows
docker run -it --workdir /data -v %cd%:/data:z immcantation/suite:4.3.0 bash

# Singularity
singularity exec -B $10x_data_2subj:/data immcantation_suite-4.3.0.sif bash
```

## Without Docker

You can run the Jupyter notebooks locally. You will need to:

* Install Jupyter and dependencies:

      * (IRkernel)[https://github.com/IRkernel/IRkernel]

      * (rpy2)[https://rpy2.github.io]

      * (bash_kernel)[https://pypi.org/project/bash_kernel]

* Install the Immcantation suite

* (Install IgBLAST)[https://changeo.readthedocs.io/en/stable/examples/igblast.html].

* Download the example data.

* Install dependencies specific for each tutorial
