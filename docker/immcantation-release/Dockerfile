FROM kleinstein/immcantation:base
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Release build of the Immcantation framework."

# Version
COPY Version.yaml /Version.yaml

# Define pipeline scripts
ENV PRESTO_ABSEQ="PrestoV5.4_AbSeq.sh" \
    CHANGEO_IGBLAST="ChangeoV3.7_IgBLASTV1.7.sh" \
    CHANGEO_CLONE="ChangeoV3.7_Clone.sh" \
    SHAZAM_THRESHOLD="ShazamV1.8_Distance.R" \
    TIGGER_GENOTYPE="TiggerV2.10_Genotype.R"

# Define R build options
ARG R_REPO="http://lib.stat.cmu.edu/R/CRAN"
ARG R_RCPP="0.12.11"
ARG R_DEPS="c('roxygen2', 'testthat', 'rmarkdown', 'knitr', 'optparse')"
ARG R_BUILD="library(devtools); install_deps(upgrade=FALSE, clean=TRUE); document(); install(build_vignettes=TRUE, clean=TRUE)"

# Setup R environment
RUN echo "options(repos='${R_REPO}', download.file.method='libcurl')" >> /usr/lib64/R/library/base/R/Rprofile \
    && Rscript -e "install.packages('devtools', clean=TRUE)" \
    && Rscript -e "devtools::install_version('Rcpp', '${R_RCPP}', quiet=FALSE, clean=TRUE)" \
    && Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"

# Install procotol data, utility scripts and pipelines
RUN PACKAGE="immcantation" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
	&& mv /tmp/${PACKAGE}/scripts/fastq2fasta.py /usr/local/bin/fastq2fasta.py \
	&& mv /tmp/${PACKAGE}/scripts/fetch_igblastdb.sh /usr/local/bin/fetch_igblastdb.sh \
	&& mv /tmp/${PACKAGE}/scripts/fetch_imgtdb.sh /usr/local/bin/fetch_imgtdb.sh \
	&& mv /tmp/${PACKAGE}/scripts/clean_imgtdb.py /usr/local/bin/clean_imgtdb.py \
	&& mv /tmp/${PACKAGE}/scripts/imgt2igblast.sh /usr/local/bin/imgt2igblast.sh \
	&& mv /tmp/${PACKAGE}/scripts/run_igblast.sh /usr/local/bin/run_igblast.sh \
	&& mv /tmp/${PACKAGE}/pipelines/${PRESTO_ABSEQ} /usr/local/bin/presto-abseq \
	&& mv /tmp/${PACKAGE}/pipelines/${CHANGEO_IGBLAST} /usr/local/bin/changeo-igblast \
	&& mv /tmp/${PACKAGE}/pipelines/${CHANGEO_CLONE} /usr/local/bin/changeo-clone \
	&& mv /tmp/${PACKAGE}/pipelines/${SHAZAM_THRESHOLD} /usr/local/bin/shazam-threshold \
	&& mv /tmp/${PACKAGE}/pipelines/${TIGGER_GENOTYPE} /usr/local/bin/tigger-genotype \
	&& mv /tmp/${PACKAGE}/protocols /usr/local/share/protocols \
	&& rm -r /tmp/${PACKAGE}

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

# Install PHYLIP
RUN PHYLIP=$(versions get -n phylip) \
    && wget -q --show-progress \
       http://evolution.gs.washington.edu/phylip/download/phylip-${PHYLIP}.tar.gz \
    && tar -zxf phylip-${PHYLIP}.tar.gz \
    && (cd phylip-${PHYLIP}/src && make -f Makefile.unx install) \
    && mv phylip-${PHYLIP}/exe/* /usr/local/bin \
    && rm -r phylip-${PHYLIP}.tar.gz phylip-${PHYLIP}

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

# Download and build references databases
RUN fetch_igblastdb.sh -o /usr/local/share/igblast \
    && fetch_imgtdb.sh -o /usr/local/share/germlines/imgt \
    && imgt2igblast.sh -i /usr/local/share/germlines/imgt -o /usr/local/share/igblast

# Setup users and permissions
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
            "  tigger-genotype\n" \
            "  shazam-threshold"
