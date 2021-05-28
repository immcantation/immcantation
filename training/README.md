# Immcantation training materials

[TOC]

## Introduction to B cell repertoire analysis 

Get a global overview of how the different tools in the Immcantation framework work together with a [Jupyter notebook](intro-lab.ipynb?viewer=nbviewer) based on the materials presented in the [webinar](https://immcantation.eventbrite.com). Use it online with 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/immcantation/immcantation-lab/master) or locally, following the instructions below. Learn more about Jupyter notebooks [here](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/). If you don’t want to execute the Jupyter notebook at this time, you can explore [the website version, with the results, here](https://kleinstein.bitbucket.io/tutorials/intro-lab/index.html).

### Learning Outcomes

* V(D)J gene annotation and novel polymorphism detection
* Clonotype assignment
* Diversity analysis
* Mutational load profiling
* Modeling of somatic hypermutation targeting
* Quantification of selection pressure

### Using the module locally with Docker

* 1.1 Pull the Immcantation Lab container image:

```
# Example: pull the development version of the lab
docker pull immcantation/lab:devel
```
    
* 1.2 Run the container:

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

* 1.3 Open your computer's internet browser and visit the url

When you visit the url from the previous step, you will start a Jupyter session in your browser.

```
# Example: http://localhost:8888/?token=18303237b2521e72f00685e4fdf754f955f82a958a8e57ec
```

* 1.4 Open the notebook

Open the notebook you want to work with. Use CTRL+Enter to execute the commands inside the cells.

For an introduction to Jupyter, visit the [official documentation site](https://jupyter-notebook.readthedocs.io/en/latest/).


### Target audience

Bioinformaticians, immunologists, biologist, scientists and students learning about immune repertoires.

### Prerequisite Skills and Knowledge Required

Familiarity with R, python and Linux.

### Domain Problem

The field of high-throughput adaptive immune receptor repertoire sequencing (AIRR-seq) has experienced significant growth in recent years, but this growth has come with considerable complexity and variety in experimental design. These complexities, combined with the high germline and somatic diversity of immunoglobulin repertoires, present analytical challenges requiring specialized methodologies. This tutorial will cover common investigative approaches and pitfalls in AIRR-seq data analysis.

### Dataset for the case study

Processed reads (input.fasta) from one healthy donor (PGP1) 3 weeks after flu vaccination (Laserson et al. (2014))

### Funding

Developed with the support of the National Library of Medicine (NIH NLM T15 LM007056) and the National Institute of Allergy and Infectious Diseases (NIH NIAID R01 AI104739).

### License

This work is licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License](https://creativecommons.org/licenses/by-nc/4.0/).


## Integration of BCR data and GEX data

Learn how to integrate BCR repertoire and gene expression analysis with this [Jupyter notebook](BCR_Seurat_tutorial.ipynb?viewer=nbviewer).