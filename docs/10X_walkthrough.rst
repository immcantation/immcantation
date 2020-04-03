10X BCR sequencing walkthrough
===========================================================================================

Overview
-------------------------------------------------------------------------------------------

This tutorial is simple walkthrough for defining B cell clonal families and building
B cell lineage trees using 10X BCR sequencing data. It is intended for users without prior experience with Immcantation. If you are familiar with Immcantation, `this page <https://changeo.readthedocs.io/en/stable/examples/10x.html>`__ may be more useful. Knowledge of basic command line usage is assumed. This tutorial will only scratch the surface of these areas. Please check out the documentation sites for the functions detailed in this tutorial before using them on your own data. For simplicity, this tutorial will use the `Immcantation Docker image <https://immcantation.readthedocs.io/en/stable/docker/intro.html>`__ which contains all necessary software. It is also possible to install the programs used separately, see `pRESTO <http://presto.readthedocs.io>`__, `Change-O <http://changeo.readthedocs.io>`__, and `Alakazam <http://alakazam.readthedocs.io>`__.


Getting started
-------------------------------------------------------------------------------------------

First, `download and unzip the example data <https://drive.google.com/open?id=17OnHtCcqV29LqyP5p8W4HR-nRnCfIirJ>`__. These data are Ig V(D)J sequences obtained PBMCs from a healthy human donor. This data was provided by 10X Genomics, and was processed using their Cell Ranger pipeline. The original data can be obtained from `here <https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.0/vdj_v1_hs_pbmc2_b?>`__.

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

Assign V(D)J gene segments
-------------------------------------------------------------------------------------------
The first steps in this pipeline are to assign V, D, and J gene segments using IgBLAST, and then convert the output into Change-O TSV format. The options of these commands differ by species. Note also that this will overwrite gene segments called by Cell Ranger. For a full treatment of this topic, see `using IgBLAST <https://changeo.readthedocs.io/en/stable/examples/igblast.html>`__ and  `parsing 10X data <https://changeo.readthedocs.io/en/stable/examples/10x.html>`__. It is also important to note that this uses the standard IMGT reference database of human alleles. To infer novel alleles and subject-specific genotypes, which would result in more accurate assignments, use `TIgGER <https://tigger.readthedocs.io/en/stable/vignettes/Tigger-Vignette/>`__.

In the Docker container, enter::

 # Assign V, D, and J segments (may take a few minutes)
 AssignGenes.py igblast -s filtered_contig.fasta -b /usr/local/share/igblast \
    --organism human --loci ig --format blast

 # Convert to AIRR format
 MakeDb.py igblast -i filtered_contig_igblast.fmt7 -s filtered_contig.fasta \
    -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_*.fasta \
    --10x filtered_contig_annotations.csv --outname 10x_igblast

This will result in the file '10x_igblast_db-pass.tab', which is the basis for the rest of the tutorial.

Identify clonal groups
-------------------------------------------------------------------------------------------
Clonal groups are B cells that descend from a common naive B cell ancestor. To group sequences into inferred clonal groups, we cluster BCR sequences that have the same heavy chain V and J segment and same junction length. We next cluster sequences with similar junction regions, using either a `defined sequence distance cutoff <https://changeo.readthedocs.io/en/stable/examples/cloning.html>`__, or an `adaptive threshold <https://scoper.readthedocs.io/en/stable/>`__. When available, we can also split clonal groups that have `differing light chain V and J segments. <https://changeo.readthedocs.io/en/stable/examples/10x.html>`__ First, separate heavy and light chain sequences::

 # Separate heavy and light chains
 ParseDb.py select -d 10x_igblast_db-pass.tab -f LOCUS -u "IGH" \
     --logic all --regex --outname heavy

 ParseDb.py select -d 10x_igblast_db-pass.tab -f LOCUS -u "IG[LK]" \
     --logic all --regex --outname light

Next, define clonal groups using a fixed heavy junction distance threshold. The threshold we use here (--dist 0.1) is not appropriate for all datasets. Check `here <https://changeo.readthedocs.io/en/stable/examples/cloning.html>`__ for instruction on determining the correct threshold for your own data. Alternatively, you could use an adaptive threshold in `SCOPer <https://scoper.readthedocs.io/en/stable/>`__. To assign clones with the distance threshold 0.1, enter::

 # Define heavy chain clones
 DefineClones.py -d heavy_parse-select.tab --act set --model ham \
     --norm len --dist 0.1 --outname heavy

Once we have defined clonal groups using heavy chains as above, we can split these groups based on whether they have differing light chain V and J gene segments::

 # Split heavy chain clones with different light chains
 light_cluster.py -d heavy_clone-pass.tab -e light_parse-select.tab \
     -o 10X_clone-pass.tab

This results in the file '10X_clone-pass.tab' which contains heavy chain sequence information derived from '10x_igblast_db-pass.tab' with an additional column 'CLONE' specifying the clonal group of the sequence.

Reconstruct heavy chain germlines
-------------------------------------------------------------------------------------------

We'll next want to reconstruct the germline heavy chain V and J segment for each clone, which will allow us to build lineage trees. Note that the options specified mask the D gene segment using 'N' nucleotides because the D gene is difficult to accurately reconstruct. For more guidance on reconstructing germlines `see this page <https://changeo.readthedocs.io/en/stable/examples/germlines.html>`__. Enter::

 # Reconstruct germline V and J sequences
 CreateGermlines.py -d 10X_clone-pass.tab -g dmask --cloned \
    -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta \
    /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHD.fasta \
    /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHJ.fasta \
    --outname 10X 

Build lineage trees
-------------------------------------------------------------------------------------------
Lineage trees represent the series of shared and unshared mutations leading from clone's germline sequence to the observed sequence data. There are multiple ways of building and visualizing these trees. Currently the simplest within Immcantation is to use `Alakazam <https://alakazam.readthedocs.io>`__, which is built around building maximum parsimony trees using PHYLIP. Alternatively, you can use `IgPhyML <https://igphyml.readthedocs.io>`__, which builds maximum likelihood trees with B cell specific models. Here, for simplicity, we use Alakazam. For more detail see Alakazam's `lineage vignette <https://alakazam.readthedocs.io/en/stable/vignettes/Lineage-Vignette/>`__

The commands in this section are meant to be entered into an R session. Open R within the Docker container using the command R. Once inside the R session, load the appropriate libraries and read in the data::

 library(alakazam)
 library(igraph)
 library(dplyr)
 
 # read data
 db <- readChangeoDb("10X_germ-pass.tab")

 # remove cells without a constant region call
 db <- filter(db, !is.na(C_CALL))

We next process clones into objects that can be used by Alakazam. This function will collapse all identical sequences within each clones, and has many options to specify which fields should be copied from the original data frame to the clone objects (i.e. text_fields)::

 # Preprocess clones
 clones <- db %>%
    group_by(CLONE) %>%
    do(CHANGEO=makeChangeoClone(., 
    text_fields=c("C_CALL", "CELL"), 
    num_fields="CONSCOUNT"))

We can now build the trees using PHYLIP. The variable 'dnapars_exec' refers to the location of the program dnapars within the Docker container::

 dnapars_exec <- "/usr/local/bin/dnapars"
 
 #build trees
 graphs <- lapply(clones$CHANGEO, buildPhylipLineage, 
      dnapars_exec=dnapars_exec, rm_temp=TRUE)

 # remove trees with < 2 sequences
 graphs[sapply(graphs, is.null)] <- NULL

Once built, we can visualize these trees using igraph. Here, we only visualize one tree, using default parameters. However, there are many ways to make more attractive lineage tree plots, detailed in Alakazam's `lineage vignette <https://alakazam.readthedocs.io/en/stable/vignettes/Lineage-Vignette/>`__. Enter into the R session::

 graph <- graphs[[1]]

 # save tree as a png image in the data directory
 png("graph.png",width=6,height=6,unit="in",res=300)
 plot(graph,layout=layout_as_tree)
 dev.off()

.. figure:: _static/graph.png
   :scale: 30 %
   :align: center
   :alt: graph

   Graph-formatted lineage tree of example clone 1.

The nodes of this tree represent observed and inferred sequences, while the edge labels represent the number of heavy chain mutations between the nodes. If you prefer more bifurcating trees, these are also detailed in Alakazam's `lineage vignette <https://alakazam.readthedocs.io/en/stable/vignettes/Lineage-Vignette/#converting-between-graph-phylo-and-newick-formats>`__.

To get the sequence attributes of the observed and inferred nodes within the tree, enter::

 attributes <- data.frame(vertex_attr(graph))