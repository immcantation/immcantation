
.. _10X-walkthrough:

10X Genomics V(D)J Sequence Analysis Tutorial
===========================================================================================

Overview
-------------------------------------------------------------------------------------------

This tutorial is basic walkthrough for defining B cell clonal families and building
B cell lineage trees using 10X BCR sequencing data. It is intended for users without prior experience with Immcantation. If you are familiar with Immcantation, `this page <https://changeo.readthedocs.io/en/stable/examples/10x.html>`__ may be more useful. Knowledge of basic command line usage is assumed. This tutorial will only scratch the surface of these areas. Please check out the documentation sites for the functions detailed in this tutorial before using them on your own data. For simplicity, this tutorial will use the `Immcantation Docker image <https://immcantation.readthedocs.io/en/stable/docker/intro.html>`__ which contains all necessary software. It is also possible to install the programs used separately, see `pRESTO <http://presto.readthedocs.io>`__, `Change-O <http://changeo.readthedocs.io>`__, and `Alakazam <http://alakazam.readthedocs.io>`__.


Getting started
-------------------------------------------------------------------------------------------

First, `download and unzip the example data <https://drive.google.com/open?id=1oRyGG5mYZBGgS7nnhjmhJpsDHjsfRz_I>`__. These data are Ig V(D)J sequences obtained PBMCs from a healthy human donor. This data is based on data provided by 10X Genomics, and processed using their Cell Ranger pipeline. We added a few sequences for instructional purposes. The original data can be obtained from `here <https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.0/vdj_v1_hs_pbmc2_b?>`__.

Second, install `Docker <https://www.docker.com/products/docker-desktop>`__ and
download the `Immcantation Docker image <https://immcantation.readthedocs.io/en/stable/docker/intro.html>`__. For some operating systems it may be necessary use super-user privledges (sudo), and/or to have 
`Docker Desktop <https://hub.docker.com/editions/community/docker-ce-desktop-windows>`__
running before entering these commands. In a terminal, enter::

 # download Immcantation Docker image (may take a few minutes)
 docker pull kleinstein/immcantation:3.0.0

Move to the directory where you've placed the example data and load it into the Docker image depending on your operating system::

 # Linux/Mac OS X, load Docker image in current directory
 docker run -it --workdir /data -v $(pwd):/data:z kleinstein/immcantation:3.0.0 bash

 # Windows, load Docker image in current directory
 docker run -it --workdir /data -v %cd%:/data:z kleinstein/immcantation:3.0.0 bash

Once inside the container, check everything is properly configured. Enter the following commands::

 BuildTrees.py --version
 # should return BuildTrees.py: 0.4.6 2019.07.19

 ls
 # should show filtered_contig_annotations.csv and filtered_contig.fasta, possibly others 

If the first command doesn't return the expected output, you probably aren't inside the right (or any) Docker container. If the second doesn't return the expected output, you may not be running the Docker image from the correct directory.

Assign gene segments and define clonal groups
-------------------------------------------------------------------------------------------

Most of the processing for 10X data can be handled by the ``changeo-10x`` script supplied in the Docker container. This script will automatically:

+ `Assign V, D, and J segments using IgBLAST <https://changeo.readthedocs.io/en/stable/examples/igblast.html>`__
+ `Convert IgBLAST output to Change-O format. <https://changeo.readthedocs.io/en/stable/examples/igblast.html#processing-the-output-of-igblast>`__
+ `Remove nonfunctional sequneces <https://changeo.readthedocs.io/en/stable/examples/filtering.html>`__. 
+ Define clonal groups (but see :ref:`next section <x-clones>`).
+ `Reconstruct heavy chain germline V and J sequences. <https://changeo.readthedocs.io/en/stable/examples/germlines.html>`__
+ Gather and compress intermediate files.

To run this script on the example dataset, enter the following command in the Docker container::

 # Run 10X processing script
 changeo-10x \
 -o . \
 -s filtered_contig.fasta \
 -a filtered_contig_annotations.csv \
 -g human -t ig -x 0.1

The options ``-s`` and ``-a`` refer to the sequence and sequence annotation file outputs from Cell Ranger, respectively. The ``-g`` option indicates species. The ``-x`` option specifies junction distance threshold used for assigning sequences into clonal clusters. For a full treatment of these individual operations, see `using IgBLAST <https://changeo.readthedocs.io/en/stable/examples/igblast.html>`__ and  `parsing 10X data <https://changeo.readthedocs.io/en/stable/examples/10x.html>`__. It is also important to note that this uses the standard `IMGT <http://www.imgt.org/>`__ reference database of human alleles. To infer novel alleles and subject-specific genotypes, which would result in more accurate assignments, see `TIgGER <https://tigger.readthedocs.io/en/stable/vignettes/Tigger-Vignette/>`__.


.. _x-clones:

Define clonal groups manually
-------------------------------------------------------------------------------------------
Clonal groups are B cells that descend from a common naive B cell ancestor. To group sequences into inferred clonal groups, we cluster BCR sequences that have the same heavy chain V and J segment and same junction length. We next cluster sequences with similar junction regions, using either a `defined sequence distance cutoff <https://changeo.readthedocs.io/en/stable/examples/cloning.html>`__, or an `adaptive threshold <https://scoper.readthedocs.io/en/stable/>`__. When available, we can also split clonal groups that have `differing light chain V and J segments. <https://changeo.readthedocs.io/en/stable/examples/10x.html>`__

In the above section, we used a predefined clonal clustering threshold of ``0.1`` using the ``-x`` option in the ``changeo-10x`` script. This is not appropriate for all datasets. Current best practice is to find the appropriate threshold for a given dataset. This can be done automatically in the ``changeo-10x`` script by specifying ``-x auto``.

Using ``-x auto`` to assign clones doesn't always work, however. If this command fails, there are other options for manually defining clones from the file ``filtered_contig_heavy_FUNCTIONAL-T.tab``. The first is by inspecting `a plot of sequence distances <https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/>`__. This is supplied in the file ``filtered_contig_threshold-plot.pdf``. You can then define clones manually using the chosen threshold (for example ``0.09``)::

 # Define heavy chain clones
 DefineClones.py -d filtered_contig_heavy_FUNCTIONAL-T.tab --act set --model ham \
     --norm len --dist 0.09 --outname filtered_contig_heavy

If the sequence distance plot is not bimodal, it may be more appropriate to instead use `SCOPer <https://scoper.readthedocs.io/en/stable/>`__ to assign clones using an adaptive threshold. Just be sure to name the output file ``filtered_contig_heavy_clone-pass.tab``.

Once we have defined clonal groups using heavy chains as above, we can split these groups based on whether they have differing light chain V and J gene segments. We also `reconstruct the heavy chain germline V and J segments <https://changeo.readthedocs.io/en/stable/examples/germlines.html>`__::

 # Split heavy chain clones with different light chains
 light_cluster.py -d filtered_contig_heavy_clone-pass.tab -e filtered_contig_light_FUNCTIONAL-T.tab \
     -o filtered_contig_heavy_clone-light.tab

 # Reconstruct germline V and J sequences
 CreateGermlines.py -d filtered_contig_heavy_clone-light.tab -g dmask --cloned \
    -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta \
    /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHD.fasta \
    /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHJ.fasta \
    --outname filtered_contig_heavy

This results in the file ``filtered_contig_heavy_germ-pass.tab`` which contains heavy chain sequence information derived from ``10x_igblast_db-pass.tab`` with an additional column ``CLONE`` specifying the clonal group of the sequence.

Build lineage trees
-------------------------------------------------------------------------------------------
Lineage trees represent the series of shared and unshared mutations leading from clone's germline sequence to the observed sequence data. There are multiple ways of building and visualizing these trees. Currently the simplest within Immcantation is to use `Alakazam <https://alakazam.readthedocs.io>`__, which is built around building maximum parsimony trees using `PHYLIP <http://evolution.genetics.washington.edu/phylip.html>`__. Alternatively, you can use `IgPhyML <https://igphyml.readthedocs.io>`__, which builds maximum likelihood trees with B cell specific models. Here, for simplicity, we use Alakazam here. For more detail see Alakazam's `lineage vignette <https://alakazam.readthedocs.io/en/stable/vignettes/Lineage-Vignette/>`__

The commands in this section are meant to be entered into an ``R`` session. Open ``R`` within the Docker container using the command ``R``. Once inside the ``R`` session, load the appropriate libraries and read in the data::

 library(alakazam)
 library(igraph)
 library(dplyr)
 
 # read data
 db <- readChangeoDb("filtered_contig_heavy_germ-pass.tab")

 # remove cells without a constant region call
 db <- filter(db, !is.na(C_CALL))

We next process clones into objects that can be used by `Alakazam <https://alakazam.readthedocs.io>`__. This function will collapse all identical sequences within each clones, and has many options to specify which fields should be copied from the original data frame to the clone objects (i.e. ``text_fields``)::

 # Preprocess clones
 clones <- db %>%
     group_by(CLONE) %>%
     do(CHANGEO=makeChangeoClone(.,
       id="CELL", text_fields=c("C_CALL"), 
         num_fields="CONSCOUNT"))

We can now build the trees using `PHYLIP <http://evolution.genetics.washington.edu/phylip.html>`__. The variable ``dnapars_exec`` refers to the location of the PHYLIP program ``dnapars`` within the Docker container::

 dnapars_exec <- "/usr/local/bin/dnapars"
 
 #build trees
 graphs <- lapply(clones$CHANGEO, buildPhylipLineage, 
      dnapars_exec=dnapars_exec, rm_temp=TRUE)

 # remove trees with < 2 sequences
 graphs[sapply(graphs, is.null)] <- NULL

Once built, we can visualize these trees using igraph. Here, we only visualize one tree, using default parameters. However, there are many ways to make more attractive lineage tree plots, detailed in Alakazam's `lineage vignette <https://alakazam.readthedocs.io/en/stable/vignettes/Lineage-Vignette/>`__. Enter into the ``R`` session::

 graph <- graphs[[1]]

 # save tree as a png image in the data directory
 png("graph.png",width=6,height=6,unit="in",res=300)
 plot(graph,layout=layout_as_tree)
 dev.off()

.. figure:: ../_static/graph.png
   :scale: 30 %
   :align: center
   :alt: graph

   Graph-formatted lineage tree of example clone 1.

The nodes of this tree represent observed and inferred sequences, while the edge labels represent the number of heavy chain mutations between the nodes. If you prefer  bifurcating trees, these are also detailed in Alakazam's `lineage vignette <https://alakazam.readthedocs.io/en/stable/vignettes/Lineage-Vignette/#converting-between-graph-phylo-and-newick-formats>`__.

To get the sequence attributes of the observed and inferred nodes within the tree, enter::

 attributes <- data.frame(vertex_attr(graph))