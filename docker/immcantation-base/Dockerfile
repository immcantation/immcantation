FROM fedora:25
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Standard environment and dependencies for the Immcantation framework."

# Bind points
VOLUME /data
VOLUME /scratch
VOLUME /software
VOLUME /oasis

# Versioning
COPY tools/versions.py /usr/local/bin/versions
COPY tools/builds.py /usr/local/bin/builds

# Update and install required packages
RUN dnf -y update && dnf install -y \
    automake \
    atlas-devel \
    biber \
    findutils \
    git \
    lapack-devel \
    libcurl-devel \
    libxml2-devel \
    mercurial \
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

# Setup environment
ENV LD_LIBRARY_PATH="/usr/lib64/R/lib:/usr/lib64/atlas"
ARG R_REPO="https://cran.cnr.berkeley.edu"
RUN echo "options(repos='${R_REPO}')" >> /usr/lib64/R/library/base/R/Rprofile