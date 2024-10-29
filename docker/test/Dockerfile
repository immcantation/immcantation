FROM immcantation/base:devel
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
              'KernSmooth', \
              'knitr', \
              'lazyeval', \
              'markdown', \
              'MASS', \
              'openxlsx', \
              'phangorn', \
              'phylotate', \
              'progress', \
              'readr', \
              'rlang', \
              'rmarkdown', \
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
RUN Rscript -e "install.packages(${R_DEPS}, clean=TRUE); devtools::install_github('javh/markr')"

# Install procotol data, utility scripts and pipelines
RUN PACKAGE="immcantation" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && mv /tmp/${PACKAGE}/scripts/* /usr/local/bin/ \
    && mv /tmp/${PACKAGE}/protocols /usr/local/share/protocols \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install muscle
RUN MUSCLE=$(versions get -n muscle) \
    && wget2 -q --force-progress --no-check-certificate \
       http://www.drive5.com/muscle/muscle${MUSCLE}_binaries.tar.gz \
	&& tar -zxf muscle${MUSCLE}_binaries.tar.gz muscle${MUSCLE}_i86linux64 \
	&& mv muscle${MUSCLE}_i86linux64 /usr/local/bin/muscle \
	&& rm muscle${MUSCLE}_binaries.tar.gz

# Install vsearch
RUN VSEARCH=$(versions get -n vsearch) \
    && wget2 -q --force-progress --no-check-certificate \
       http://github.com/torognes/vsearch/releases/download/v${VSEARCH}/vsearch-${VSEARCH}-linux-x86_64.tar.gz \
	&& tar -zxf vsearch-${VSEARCH}-linux-x86_64.tar.gz \
	&& mv vsearch-${VSEARCH}-linux-x86_64/bin/vsearch /usr/local/bin/vsearch \
	&& ln -s /usr/local/bin/vsearch /usr/local/bin/usearch \
	&& rm -r vsearch-${VSEARCH}-linux-x86_64.tar.gz vsearch-${VSEARCH}-linux-x86_64

# Install BLAST+ executables
RUN BLAST=$(versions get -n blast) \
    && wget2 -q --force-progress --no-check-certificate \
       https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST}/ncbi-blast-${BLAST}+-x64-linux.tar.gz \
    && tar -zxf ncbi-blast-${BLAST}+-x64-linux.tar.gz \
    && mv ncbi-blast-${BLAST}+/bin/* /usr/local/bin \
    && rm -r ncbi-blast-${BLAST}+-x64-linux.tar.gz ncbi-blast-${BLAST}+

# Install IgBLAST executables
RUN IGBLAST=$(versions get -n igblast) \
    && wget2 -q --force-progress --no-check-certificate \
       https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/${IGBLAST}/ncbi-igblast-${IGBLAST}-x64-linux.tar.gz \
    && tar -zxf ncbi-igblast-${IGBLAST}-x64-linux.tar.gz \
    && mv ncbi-igblast-${IGBLAST}/bin/* /usr/local/bin \
    && mkdir -p /usr/local/share/igblast \
    && fetch_igblastdb.sh -o /usr/local/share/igblast \
    && mv ncbi-igblast-${IGBLAST}/internal_data ncbi-igblast-${IGBLAST}/optional_file /usr/local/share/igblast \
    && rm -r ncbi-igblast-${IGBLAST}-x64-linux.tar.gz ncbi-igblast-${IGBLAST}

# Install CD-HIT
RUN CDHIT=$(versions get -n cd-hit) \
    && wget2 -q --force-progress --no-check-certificate \
       https://github.com/weizhongli/cdhit/archive/V${CDHIT}.tar.gz \
    && tar -zxf V${CDHIT}.tar.gz \
    && (cd cdhit-${CDHIT} && make && make install) \
    && rm -r V${CDHIT}.tar.gz cdhit-${CDHIT}

# Install PHYLIP
# Temporary fix for makebox, added following the patch for debian:
# https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=1066595
RUN PHYLIP=$(versions get -n phylip) \
    && wget2 -q --force-progress --no-check-certificate \
       https://phylipweb.github.io/phylip/download/phylip-${PHYLIP}.tar.gz \
    && tar -zxf phylip-${PHYLIP}.tar.gz \
    && (cd phylip-${PHYLIP}/src && sed -i '/double computeAngle/a void makebox(char *, double *, double *, double *, long);' draw.h && make -f Makefile.unx install CFLAGS='-fcommon') \
    && mv phylip-${PHYLIP}/exe/* /usr/local/bin \
    && rm -r phylip-${PHYLIP}.tar.gz phylip-${PHYLIP}

# Install IgPhyML
RUN IGPHYML=$(versions get -n igphyml) \
    && wget2 -q --force-progress --no-check-certificate \
       https://bitbucket.org/kleinstein/igphyml/get/${IGPHYML}.tar.gz -O igphyml-${IGPHYML}.tar.gz \
    && mkdir -p /usr/local/share/igphyml \
    && tar -zxf igphyml-${IGPHYML}.tar.gz -C /usr/local/share/igphyml --strip-components 1 \
    && (cd /usr/local/share/igphyml && ./make_phyml_blas_omp) \
    && rm -r igphyml-${IGPHYML}.tar.gz

# Install AIRR reference libraries
RUN AIRR_PY=$(versions get -n airr-py) \
    && AIRR_R=$(versions get -n airr-r) \
    && pip3 install airr==${AIRR_PY} \
    && Rscript -e "devtools::install_version('airr', '${AIRR_R}', quiet=FALSE, clean=TRUE)"

# Download and build references databases
RUN fetch_imgtdb.sh -o /usr/local/share/germlines/imgt \
    && imgt2igblast.sh -i /usr/local/share/germlines/imgt -o /usr/local/share/igblast \
    && fetch_phix.sh -o /usr/local/share/phix

# Patch IgBLAST database
#RUN wget2 -q --force-progress --no-check-certificate \
#    ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/patch/optional_file/mouse_gl.aux \
#    && mv -f mouse_gl.aux /usr/local/share/igblast/optional_file

# Set commands
CMD echo -e "Report version information:\n" \
            "  versions report\n" \
            "Report build details:\n" \
            "  builds report\n"
