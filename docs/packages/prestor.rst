.. _prestoR:

prestoR
================================================================================

The presto report package (prestoR) is an R package for generating
quality control plots from `pRESTO <http://presto.readthedocs.io>`_ log tables.

:download:`Example Report <../_static/example_report.pdf>`

Download & Installation
--------------------------------------------------------------------------------

`prestor` is current not available from CRAN and must be installed from the
GitHub repo directly by first cloning the repository:

`https://github.com/immcantation/prestor <https://github.com/immcantation/prestor>`_

Then build using the following R commands from the package root::

    install.packages(c("devtools", "roxygen2"))
    library(devtools)
    install_deps(dependencies=T)
    document()
    install()

Alternatively, you can install directly form the GitHub repository, but this
will not build the documentation::

    library(devtools)
    install_github("immcantation/prestor@master")

Documentation
--------------------------------------------------------------------------------

For an index of available functions see::

    help(package="prestor")

For some common tasks, see the following help pages:

====================  ===========================================================
Function              Description
====================  ===========================================================
buildReport           Generate a presto pipeline report
loadConsoleLog	      Parse console output from a pRESTO pipeline
loadLogTable	      Parse tabled log output from pRESTO tools
pdfReport	          R Markdown to PDF format for pRESTO reports
plotAlignSets	      Plot AlignSets log table
plotAssemblePairs	  Plot AssemblePairs log table
plotBuildConsensus	  Plot BuildConsensus log table
plotConsoleLog	      Plot console output from a pRESTO pipeline
plotFilterSeq	      Plot FilterSeq log table
plotMaskPrimers	      Plot MaskPrimer log table
plotParseHeaders	  Plot ParseHeaders log table
report_abseq3         Generate a report for an AbSeq V3 pRESTO pipeline script
====================  ===========================================================