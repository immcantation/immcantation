FROM kleinstein/immcantation:base
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Test environment for the Immcantation framework."

# Version
COPY Version.yaml /Version.yaml
RUN builds write -n date -v "$(date +'%Y-%m-%d %T %Z')"

# Define R package installation parameters
# ARG R_RCPP="0.12.16"
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
              'gtools', \
              'hexbin', \
              'igraph', \
              'iterators', \
              'kedd', \
              'KernSmooth', \
              'knitr', \
              'lazyeval', \
              'MASS', \
              'openxlsx', \
              'progress', \
              'readr', \
              'rlang', \
              'rversions', \
              'scales', \
              'SDMTools', \
              'seqinr', \
              'stringi', \
              'tidyr', \
              'versions', \
              'yaml')"

# Setup R environment
#RUN Rscript -e "install.packages('devtools', clean=TRUE)" \
#    && Rscript -e "devtools::install_version('Rcpp', '${R_RCPP}', quiet=FALSE, clean=TRUE)" \
#    && Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"
RUN Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"

# Install muscle
RUN MUSCLE=$(versions get -n muscle) \
    && wget -q --show-progress --no-check-certificate \
       http://www.drive5.com/muscle/muscle${MUSCLE}_binaries.tar.gz \
	&& tar -zxf muscle${MUSCLE}_binaries.tar.gz muscle${MUSCLE}_i86linux64 \
	&& mv muscle${MUSCLE}_i86linux64 /usr/local/bin/muscle \
	&& rm muscle${MUSCLE}_binaries.tar.gz

# Install vsearch
RUN VSEARCH=$(versions get -n vsearch) \
    && wget -q --show-progress --no-check-certificate \
       http://github.com/torognes/vsearch/releases/download/v${VSEARCH}/vsearch-${VSEARCH}-linux-x86_64.tar.gz \
	&& tar -zxf vsearch-${VSEARCH}-linux-x86_64.tar.gz \
	&& mv vsearch-${VSEARCH}-linux-x86_64/bin/vsearch /usr/local/bin/vsearch \
	&& ln -s /usr/local/bin/vsearch /usr/local/bin/usearch \
	&& rm -r vsearch-${VSEARCH}-linux-x86_64.tar.gz vsearch-${VSEARCH}-linux-x86_64

# Install BLAST+ executables
RUN BLAST=$(versions get -n blast) \
    && wget -q --show-progress --no-check-certificate \
       ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST}/ncbi-blast-${BLAST}+-x64-linux.tar.gz \
    && tar -zxf ncbi-blast-${BLAST}+-x64-linux.tar.gz \
    && mv ncbi-blast-${BLAST}+/bin/* /usr/local/bin \
    && rm -r ncbi-blast-${BLAST}+-x64-linux.tar.gz ncbi-blast-${BLAST}+

# Install IgBLAST executables
RUN IGBLAST=$(versions get -n igblast) \
    && wget -q --show-progress --no-check-certificate \
       ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST}/ncbi-igblast-${IGBLAST}-x64-linux.tar.gz \
    && tar -zxf ncbi-igblast-${IGBLAST}-x64-linux.tar.gz \
    && mv ncbi-igblast-${IGBLAST}/bin/* /usr/local/bin \
    && rm -r ncbi-igblast-${IGBLAST}-x64-linux.tar.gz ncbi-igblast-${IGBLAST}

# Install CD-HIT
RUN CDHIT=$(versions get -n cd-hit) \
    && wget -q --show-progress --no-check-certificate \
       https://github.com/weizhongli/cdhit/archive/V${CDHIT}.tar.gz \
    && tar -zxf V${CDHIT}.tar.gz \
    && (cd cdhit-${CDHIT} && make && make install) \
    && rm -r V${CDHIT}.tar.gz cdhit-${CDHIT}

# Install PHYLIP
RUN PHYLIP=$(versions get -n phylip) \
    && wget -q --show-progress --no-check-certificate \
       http://evolution.gs.washington.edu/phylip/download/phylip-${PHYLIP}.tar.gz \
    && tar -zxf phylip-${PHYLIP}.tar.gz \
    && (cd phylip-${PHYLIP}/src && make -f Makefile.unx install) \
    && mv phylip-${PHYLIP}/exe/* /usr/local/bin \
    && rm -r phylip-${PHYLIP}.tar.gz phylip-${PHYLIP}

# Install IgPhyML
RUN IGPHYML=$(versions get -n igphyml) \
    && wget -q --show-progress --no-check-certificate \
       https://bitbucket.org/kbhoehn/igphyml/get/${IGPHYML}.tar.gz -O igphyml-${IGPHYML}.tar.gz \
    && mkdir -p /usr/local/share/igphyml \
    && tar -zxf igphyml-${IGPHYML}.tar.gz -C /usr/local/share/igphyml --strip-components 1 \
    && (cd /usr/local/share/igphyml && ./make_phyml_blas_omp) \
    && rm -r igphyml-${IGPHYML}.tar.gz

# Set commands
CMD echo -e "Report version information:\n" \
            "  versions report\n" \
            "Report build details:\n" \
            "  builds report\n"
