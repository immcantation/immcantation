FROM fedora:42
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
    fontconfig-devel \
    fribidi-devel \
    git \
    glpk-devel \
    gmp-devel \
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
    procps \
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
    R-devtools \
    redland-devel \    
    sudo \
    tar \
    texlive-pdfcrop \
    texlive-titling \
    v8-devel \
    wget \
    wget2 \
    && dnf clean all

# System environment variables
ENV LD_LIBRARY_PATH="/usr/lib64/R/lib:/usr/lib64/atlas"
#ENV R_LIBS_SITE="/usr/lib64/R/library:/usr/share/R/library"
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 20

# R build setup
ARG R_REPO="https://cloud.r-project.org/"
ARG R_DEPS_REINSTALL="c('callr', \
                        'crayon', \
                        'desc', \
                        'lifecycle', \
                        'pkgbuild', \
                        'prettyunits', \
                        'R6', \
                        'rcmdcheck', \
                        'rprojroot', \
                        'sessioninfo', \
                        'withr', \
                        'xopen')"
ARG R_DEPS="c('BiocManager', \
              'httpuv', \
              'httr', \
              'knitr', \
              'Matrix', \
              'memoise', \
              'optparse', \
              'png', \
              'ragg', \
              'rmarkdown', \
              'roxygen2', \
              'stringr', \
              'testthat', \
              'V8', \
              'withr')"
ARG R_BIOC_DEPS="c('Biostrings', \
              'ComplexHeatmap', \
              'DECIPHER', \
              'GenomicAlignments', \
              'ggmsa', \
              'ggtree', \
              'IRanges')"
# Install R packages
# - Force install of Rcpp from CRAN first to avoid dependency issues
# - Temporary fix due to compatibility issues with enchantr:
#   - Lock ggplot2 version 3.5.2. Set dependency to only 'Depends', 'Imports', and 'LinkingTo' to avoid upgrading ggplot2 itself (?!)
#     It seems that some dependencies of ggplot2 suggest to upgrade ggplot2.
#   - Need ggiraph version to 0.9.0 to be compatible with ggplot2 3.5.2
#   - Using Bioconductor 3.21 (3.22 upgrades ggplot2)
RUN echo "options(repos='${R_REPO}')" >> /usr/lib64/R/library/base/R/Rprofile \
    && mkdir -p /usr/share/doc/R/html
RUN Rscript -e "for (lib in .libPaths()) { \
                    installed_in_lib <- rownames(installed.packages(lib.loc=lib));\
                    to_remove <- intersect(${R_DEPS_REINSTALL}, installed_in_lib);\
                    if (length(to_remove) > 0) { remove.packages(to_remove, lib=lib) }\
                }"
RUN Rscript -e "install.packages(${R_DEPS_REINSTALL}, clean=TRUE)"
RUN Rscript -e "devtools::install_cran('Rcpp', upgrade='never', dependencies=T, force=T)"
RUN Rscript -e "devtools::install_version('ggplot2', upgrade='never', dependencies=c('Depends', 'Imports', 'LinkingTo'), version='3.5.2');library('ggplot2')"
RUN Rscript -e "devtools::install_version('ggiraph', version='0.9.0', upgrade='never');library('ggiraph')"
RUN Rscript -e "devtools::install_cran(${R_DEPS}, upgrade='never', dependencies=c('Depends', 'Imports', 'LinkingTo'))"
RUN Rscript -e " \
    BiocManager::install(version='3.21', ask=FALSE, update=FALSE); \
    BiocManager::install(${R_BIOC_DEPS}, update=FALSE, ask=FALSE, force=FALSE); \
    for (p in c(${R_DEPS}, ${R_BIOC_DEPS})) {library(p,character.only=T)}; \
    "
RUN touch /.Rprofile \
    && touch /.Renviron \
    && Rscript -e "cat(paste0('R_LIBS_SITE=\'',paste(.libPaths(),collapse=':'),'\''), sep='\n')" >> /.Renviron