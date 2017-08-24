FROM kleinstein/immcantation:base
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Test environment for the Immcantation framework."

# Define R package installation parameters
ARG R_REPO="http://lib.stat.cmu.edu/R/CRAN"
ARG R_RCPP="0.12.12"
ARG R_DEPS="c('ape', \
              'bibtex', \
              'captioner', \
              'data.table', \
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
RUN echo "options(repos='${R_REPO}', download.file.method='libcurl')" >> /usr/lib64/R/library/base/R/Rprofile \
    && Rscript -e "install.packages('devtools', clean=TRUE)" \
    && Rscript -e "devtools::install_version('Rcpp', '${R_RCPP}', quiet=FALSE, clean=TRUE)" \
    && Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"