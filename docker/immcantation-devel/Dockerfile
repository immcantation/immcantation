FROM javh/immcantation:base

# Labels
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Development build of the Immcantation framework."

# Sources
ENV IMMCANT_TAG="tip"

# Pipeline scripts
ENV PRESTO_ABSEQ="PrestoV5.3_AbSeqV3.sh" \
    CHANGEO_IGBLAST="ChangeoV3.4_IgBLASTV1.6.sh" \
    CHANGEO_CLONE="ChangeoV3.4_Clone.sh" \
    SHAZAM_THRESHOLD="ShazamV1.7_Distance.R"

# Define R package dependencies
ENV R_DEPS "c('ape', \
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

# Install current default branch of presto and changeo
RUN pip3 install --no-cache-dir hg+https://bitbucket.org/kleinstein/presto#default\
    && pip3 install --no-cache-dir hg+https://bitbucket.org/kleinstein/changeo#default

# Set default CRAN repo
RUN echo "options(repos=c(CRAN='http://lib.stat.cmu.edu/R/CRAN/'), download.file.method='libcurl')" \
	>> /usr/lib64/R/library/base/R/Rprofile

# Install R package dependencies
RUN Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"

# Install current default branch of alakazam, shazam, tigger and prestor
RUN Rscript -e "library(devtools); install_bitbucket('kleinstein/alakazam@default')" \
    && Rscript -e "library(devtools); install_bitbucket('kleinstein/shazam@default')" \
    && Rscript -e "library(devtools); install_bitbucket('kleinstein/tigger@default')" \
    && Rscript -e "library(devtools); install_bitbucket('javh/prototype-prestor@default')"

# Install pipeline scripts
#COPY pipelines/${PRESTO_ABSEQ} /usr/local/bin/presto-abseq
#COPY pipelines/${CHANGEO_IGBLAST} /usr/local/bin/changeo-igblast
#COPY pipelines/${CHANGEO_CLONE} /usr/local/bin/changeo-clone

RUN wget -q --show-progress \
    https://bitbucket.org/kleinstein/immcantation/get/${IMMCANT_TAG}.tar.gz -O immcantation.tar.gz \
	&& mkdir -p immcantation; tar --wildcards --strip-components=1 -C immcantation -zxf immcantation.tar.gz \
	&& mv immcantation/pipelines/${PRESTO_ABSEQ} /usr/local/bin/presto-abseq \
	&& mv immcantation/pipelines/${CHANGEO_IGBLAST} /usr/local/bin/changeo-igblast \
	&& mv immcantation/pipelines/${CHANGEO_CLONE} /usr/local/bin/changeo-clone \
	&& mv immcantation/pipelines/${SHAZAM_THRESHOLD} /usr/local/bin/shazam-threshold \
	&& rm -r immcantation immcantation.tar.gz

# Setup users and permissions
VOLUME /data
VOLUME /scratch
VOLUME /software
RUN useradd -u 1000 -g users magus \
    && echo "umask 000" >> /home/magus/.bashrc
USER magus

# Set commands
CMD echo -e "The following pipelines are available:\n" \
            "  presto-abseq\n" \
            "  changeo-igblast\n" \
            "  changeo-clone\n" \
            "  shazam-threshold"
