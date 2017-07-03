FROM javh/immcantation:base
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Release build of the Immcantation framework."

# Define pipeline scripts
ENV PRESTO_ABSEQ="PrestoV5.3_AbSeqV3.sh" \
    CHANGEO_IGBLAST="ChangeoV3.4_IgBLASTV1.6.sh" \
    CHANGEO_CLONE="ChangeoV3.4_Clone.sh" \
    SHAZAM_THRESHOLD="ShazamV1.7_Distance.R"

# Setup R build environment
ENV R_REPO="c(CRAN='http://lib.stat.cmu.edu/R/CRAN')" \
    R_RCPP="0.12.10" \
    R_DEPS="c('roxygen2', 'testthat', 'rmarkdown', 'knitr', 'optparse')" \
    R_BUILD="library(devtools); install_deps(upgrade=FALSE, clean=TRUE); document(); install(build_vignettes=TRUE, clean=TRUE)"

RUN echo "options(repos=${R_REPO}, download.file.method='libcurl')" >> /usr/lib64/R/library/base/R/Rprofile \
    && Rscript -e "install.packages('devtools', clean=TRUE)" \
    && Rscript -e "devtools::install_version('Rcpp', '${R_RCPP}', quiet=FALSE, clean=TRUE)" \
    && Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"

# Install pipeline scripts
RUN PACKAGE="immcantation" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
	&& mv /tmp/${PACKAGE}/pipelines/${PRESTO_ABSEQ} /usr/local/bin/presto-abseq \
	&& mv /tmp/${PACKAGE}/pipelines/${CHANGEO_IGBLAST} /usr/local/bin/changeo-igblast \
	&& mv /tmp/${PACKAGE}/pipelines/${CHANGEO_CLONE} /usr/local/bin/changeo-clone \
	&& mv /tmp/${PACKAGE}/pipelines/${SHAZAM_THRESHOLD} /usr/local/bin/shazam-threshold \
	&& rm -r /tmp/${PACKAGE}

# Install presto
RUN PACKAGE="presto" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
	&& pip3 install --no-cache-dir /tmp/${PACKAGE} \
	&& rm -r /tmp/${PACKAGE}

# Install changeo
RUN PACKAGE="changeo" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
	&& pip3 install --no-cache-dir /tmp/${PACKAGE} \
	&& rm -r /tmp/${PACKAGE}

# Install alakazam
RUN PACKAGE="alakazam" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && Rscript -e "setwd('/tmp/${PACKAGE}'); ${R_BUILD}"

# Install shazam
RUN PACKAGE="shazam" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && Rscript -e "setwd('/tmp/${PACKAGE}'); ${R_BUILD}" \
    && rm -r /tmp/${PACKAGE}

# Install tigger
RUN PACKAGE="tigger" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && Rscript -e "setwd('/tmp/${PACKAGE}'); ${R_BUILD}" \
    && rm -r /tmp/${PACKAGE}

# Install prestor
RUN PACKAGE="prestor" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone https://bitbucket.org/javh/prototype-prestor /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && Rscript -e "setwd('/tmp/${PACKAGE}'); ${R_BUILD}" \
    && rm -r /tmp/${PACKAGE}

# Setup users and permissions
VOLUME /data
VOLUME /scratch
VOLUME /software
RUN useradd -u 1000 -g users magus \
    && echo "umask 000" >> /home/magus/.bashrc
USER magus

# Set commands
CMD echo -e "Report version information:\n" \
            "  versions report\n" \
            "Available pipeline commands:\n" \
            "  presto-abseq\n" \
            "  changeo-igblast\n" \
            "  changeo-clone\n" \
            "  shazam-threshold"
