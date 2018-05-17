FROM kleinstein/immcantation:base
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Test environment for the Immcantation framework."

# Version
COPY Version.yaml /Version.yaml
RUN builds write -n date -v "$(date +'%Y-%m-%d %T %Z')"

# Define R package installation parameters
ARG R_REPO="http://lib.stat.cmu.edu/R/CRAN"
ARG R_RCPP="0.12.16"
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
              'gridExtra', \
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

# Setup R environment
RUN echo "options(repos='${R_REPO}', download.file.method='libcurl')" >> /usr/lib64/R/library/base/R/Rprofile \
    && Rscript -e "install.packages('devtools', clean=TRUE)" \
    && Rscript -e "devtools::install_version('Rcpp', '${R_RCPP}', quiet=FALSE, clean=TRUE)" \
    && Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"

# Install muscle
RUN MUSCLE=$(versions get -n muscle) \
    && wget -q --show-progress \
       http://www.drive5.com/muscle/downloads${MUSCLE}/muscle${MUSCLE}_i86linux64.tar.gz \
	&& tar -zxf muscle${MUSCLE}_i86linux64.tar.gz \
	&& mv muscle${MUSCLE}_i86linux64 /usr/local/bin/muscle \
	&& rm muscle${MUSCLE}_i86linux64.tar.gz

# Install vsearch
RUN VSEARCH=$(versions get -n vsearch) \
    && wget -q --show-progress \
       http://github.com/torognes/vsearch/releases/download/v${VSEARCH}/vsearch-${VSEARCH}-linux-x86_64.tar.gz \
	&& tar -zxf vsearch-${VSEARCH}-linux-x86_64.tar.gz \
	&& mv vsearch-${VSEARCH}-linux-x86_64/bin/vsearch /usr/local/bin/vsearch \
	&& ln -s /usr/local/bin/vsearch /usr/local/bin/usearch \
	&& rm -r vsearch-${VSEARCH}-linux-x86_64.tar.gz vsearch-${VSEARCH}-linux-x86_64

# Install BLAST+ executables
RUN BLAST=$(versions get -n blast) \
    && wget -q --show-progress \
       ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST}/ncbi-blast-${BLAST}+-x64-linux.tar.gz \
    && tar -zxf ncbi-blast-${BLAST}+-x64-linux.tar.gz \
    && mv ncbi-blast-${BLAST}+/bin/* /usr/local/bin \
    && rm -r ncbi-blast-${BLAST}+-x64-linux.tar.gz ncbi-blast-${BLAST}+

# Install IgBLAST executables
RUN IGBLAST=$(versions get -n igblast) \
    && wget -q --show-progress \
       ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST}/ncbi-igblast-${IGBLAST}-x64-linux.tar.gz \
    && tar -zxf ncbi-igblast-${IGBLAST}-x64-linux.tar.gz \
    && mv ncbi-igblast-${IGBLAST}/bin/* /usr/local/bin \
    && rm -r ncbi-igblast-${IGBLAST}-x64-linux.tar.gz ncbi-igblast-${IGBLAST}

# Install CD-HIT
RUN CDHIT=$(versions get -n cd-hit) \
    && wget -q --show-progress \
       https://github.com/weizhongli/cdhit/archive/V${CDHIT}.tar.gz \
    && tar -zxf V${CDHIT}.tar.gz \
    && (cd cdhit-${CDHIT} && make && make install) \
    && rm -r V${CDHIT}.tar.gz cdhit-${CDHIT}

# Install PHYLIP
RUN PHYLIP=$(versions get -n phylip) \
    && wget -q --show-progress \
       http://evolution.gs.washington.edu/phylip/download/phylip-${PHYLIP}.tar.gz \
    && tar -zxf phylip-${PHYLIP}.tar.gz \
    && (cd phylip-${PHYLIP}/src && make -f Makefile.unx install) \
    && mv phylip-${PHYLIP}/exe/* /usr/local/bin \
    && rm -r phylip-${PHYLIP}.tar.gz phylip-${PHYLIP}

# Install IgPhyML
RUN PACKAGE="igphyml" \
    && rm -rf /usr/local/share/${PACKAGE} \
    && git clone https://bitbucket.org/kbhoehn/${PACKAGE} /usr/local/share/${PACKAGE} \
	&& (cd /usr/local/share/${PACKAGE} && ./make_phyml_blas_omp_fedora) \
	&& builds write -n ${PACKAGE} -v $(git -C /usr/local/share/${PACKAGE} show --pretty=%h | head -n 1)

# Set commands
CMD echo -e "Report version information:\n" \
            "  versions report\n" \
            "Report build details:\n" \
            "  builds report\n"
