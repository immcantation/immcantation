FROM fedora:29
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
    lapack-devel \
    libcurl-devel \
    libidn \
    libnsl \
    libxml2-devel \
    mercurial \
    mesa-libGL-devel \
    mesa-libGLU-devel \
    openssl-devel \
    pandoc \
    pandoc-citeproc \
    perl \
    python3 \
    python3-biopython \
    python3-hglib \
    python3-numpy \
    python3-pandas \
    python3-PyYAML \
    python3-setuptools \
    python3-scipy \
    R-core \
    R-core-devel \
    sudo \
    tar \
    texlive-scheme-basic \
    texlive-biblatex \
    texlive-titling \
    wget \
    && dnf clean all

# System environment variables
ENV LD_LIBRARY_PATH="/usr/lib64/R/lib:/usr/lib64/atlas"

# R build setup
ARG R_REPO="http://cran.cnr.berkeley.edu"
ARG R_DEPS="c('Rcpp', \
              'devtools', \
              'roxygen2', \
              'testthat', \
              'pkgbuild', \
              'rmarkdown', \
              'knitr', \
              'optparse')"
RUN echo "options(repos='${R_REPO}')" >> /usr/lib64/R/library/base/R/Rprofile \
    && mkdir -p /usr/share/doc/R/html \
	&& cp /usr/lib64/R/library/stats/html/R.css /usr/share/doc/R/html/R.css \
	&& Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"
