FROM fedora:38
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Standard environment and dependencies for the Immcantation framework."

# Bind points
VOLUME /data
VOLUME /scratch
VOLUME /software
VOLUME /oasis

# Tools
COPY tools/versions.py /usr/local/bin/versions
COPY tools/builds.py /usr/local/bin/builds
COPY tools/pipelines.py /usr/local/bin/pipelines
COPY tools/rinstall.R /usr/local/bin/rinstall

# Disable https on repos
RUN (cd /etc/yum.repos.d; sed -i 's/https/http/g' *.repo)

# Update and install required packages
RUN dnf -y update && dnf install -y \
    automake \
    atlas-devel \
    biber \
    findutils \
    git \
    langpacks-en \
    lapack-devel \
    libcurl-devel \
    libgit2-devel \
    libidn \
    libidn1.34 \
    libnsl \
    libsodium-devel \
    libssh2 \
    libuv \
    libxml2-devel \
    mercurial \
    mesa-libGL-devel \
    mesa-libGLU-devel \
    openssl-devel \
    pandoc \
    perl \
    proj-devel \    
    python3 \
    python3-biopython \
    python3-GitPython \
    python3-hglib \
    python3-numpy \
    python3-pandas \
    python3-pip \
    python3-PyYAML \
    python3-setuptools \
    python3-scipy \
    R-core \
    R-core-devel \
    R-rgeos \
    R-ragg \
    R-textshaping \
    redland-devel \    
    sudo \
    tar \
    texlive-titling \
    wget \
    && dnf clean all

# System environment variables
ENV LD_LIBRARY_PATH="/usr/lib64/R/lib:/usr/lib64/atlas"
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 20

# R build setup
ARG R_REPO="http://lib.stat.cmu.edu/R/CRAN/"
ARG R_DEPS="c('BiocManager', \
              'devtools', \
              'knitr', \
              'Matrix', \
              'optparse', \
              'pkgbuild', \
              'Rcpp', \
              'rmarkdown', \
              'roxygen2', \
              'testthat')"
ARG R_BIOC_DEPS="c('Biostrings', \
              'DECIPHER', \
              'GenomicAlignments', \
              'ggmsa', \
              'ggtree', \
              'IRanges')"           
RUN echo "options(repos='${R_REPO}')" >> /usr/lib64/R/library/base/R/Rprofile \
    && mkdir -p /usr/share/doc/R/html \
    && Rscript -e "install.packages(${R_DEPS}, clean=TRUE)" \
    && Rscript -e "BiocManager::install(${R_BIOC_DEPS},update=F,ask=F,method='wget')"
RUN touch /.Rprofile \
    && touch /.Renviron

