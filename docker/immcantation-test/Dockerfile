FROM kleinstein/immcantation:base
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Test environment for the Immcantation framework."

# Define R package installation parameters
ENV R_REPO="c(CRAN='http://lib.stat.cmu.edu/R/CRAN')" \
    R_DEPS="c('ape', \
              'bibtex', \
              'captioner', \
              'data.table', \
              'devtools', \
              'diagram', \
              'diptest', \
              'doParallel', \
              'dplyr', \
              'foreach', \
              'ggplot2', \
              'hexbin', \
              'igraph', \
              'iterators', \
              'kedd', \
              'KernSmooth', \
              'knitr', \
              'lazyeval', \
              'MASS', \
              'openxlsx', \
              'optparse', \
              'progress', \
              'Rcpp', \
              'readr', \
              'rmarkdown', \
              'roxygen2', \
              'rversions', \
              'scales', \
              'SDMTools', \
              'seqinr', \
              'stringi', \
              'testthat', \
              'tidyr', \
              'versions', \
              'yaml')"

# Install R package dependencies
RUN echo "options(repos=${R_REPO}, download.file.method='libcurl')" >> /usr/lib64/R/library/base/R/Rprofile \
    && Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"
