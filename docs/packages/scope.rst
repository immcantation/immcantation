SCOPe
================================================================================

SCOPe (Spectral Clustering for clOne Partitioning) provides a computational
framework for unsupervised identification B cell clones from adaptive immune
receptor repertoire sequencing (AIRR-Seq) datasets.

This method performs spectral clustering of the B cell receptor (BCR) junction
region within groups of BCR sequences sharing the same V gene, J gene, and
junction length. Rather than a fixed threshold, SCOPe uses an adaptive threshold
for clustering sequences to determine the local sequence neighborhood, which
offers an improvement in both the sensitivity and specificity over a simple
fixed threshold for all junction lengths.

Download & Installation
--------------------------------------------------------------------------------

`scope` is current not available from CRAN and must be installed from the
bitbucket repo directly by first cloning the bitbucket repository:

`http://bitbucket.org/kleinstein/scope <https://bitbucket.org/kleinstein/scope>`_

Then build using the following R commands from the package root::

    install.packages(c("devtools", "roxygen2"))
    library(devtools)
    install_deps(dependencies=T)
    document()
    install()

Alternatively, you can install directly form the bitbucket repository, but this
will not build the documentation::

    library(devtools)
    install_bitbucket("kleinstein/scope@default")

Documentation
--------------------------------------------------------------------------------

For clustering sequences into clonal groups see::

    help("defineClonesScope", package="scope")

For summary statistics an visualization of the clonal clustering results see::

    help("clonesAnalysis", package="scope")

How to Cite
--------------------------------------------------------------------------------

    **Nouri N and Kleinstein SH**. A spectral clustering-based method for
    identifying clones from high-throughput B cell repertoire sequencing data.
    *Bioinformatics, (in press).*
