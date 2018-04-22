prestoR
================================================================================

The presto report package (prestoR) is an R package for generating
quality control plots from `pRESTO <http://presto.readthedocs.io>`_ log tables.

:download:`Example Report <../_static/example_report.pdf>`

Download & Installation
--------------------------------------------------------------------------------

`prestor` is current not available from CRAN and must be installed from the
bitbucket repo directly by first cloning the bitbucket repository:

`http://bitbucket.org/javh/prototype-prestor <https://bitbucket.org/javh/prototype-prestor>`_

Then build using the following R commands from the package root::

    install.packages(c("devtools", "roxygen2"))
    library(devtools)
    install_deps(dependencies=T)
    document()
    install()

Alternatively, you can install directly form the bitbucket repository, but this
will not build the documentation::

    library(devtools)
    install_bitbucket("javh/prototype-prestor@default")

Documentation
--------------------------------------------------------------------------------

For an index of available functions see::

    help(package="prestor")

For some common tasks see the following help pages:

====================  ===========================================================
Function              Description
====================  ===========================================================
report_abseq3         Generate a report for an AbSeq V3 pRESTO pipeline script
loadConsoleLog	      Parse console output from a pRESTO pipeline
loadLogTable	      Parse tabled log output from pRESTO tools
plotConsoleLog	      Plot console output from a pRESTO pipeline
plotAssemblePairs	  Plot AssemblePairs log table
plotBuildConsensus	  Plot BuildConsensus log table
plotFilterSeq	      Plot FilterSeq log table
plotMaskPrimers	      Plot MaskPrimer log table
plotParseHeaders	  Plot ParseHeaders log table
====================  ===========================================================
