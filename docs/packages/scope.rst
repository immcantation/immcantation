SCOPe
================================================================================

Provides a computational framework for unsupervised identification B cell
clones from adaptive immune receptor repertoire sequencing (AIRR-Seq) datasets.
This method is based on spectral clustering of the junction sequences of B cell
receptors (BCRs, Immunoglobulins) that share the same V gene, J gene and
junction length.

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
