FROM kleinstein/immcantation:base
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Release build of the Immcantation framework."

# Version
COPY Version.yaml /Version.yaml
RUN builds write -n date -v "$(date +'%Y-%m-%d %T %Z')"

# Define pipeline scripts
ARG PRESTO_ABSEQ="pipelines/presto-abseq.sh"
ARG CHANGEO_IGBLAST="pipelines/changeo-igblast.sh"
ARG CHANGEO_CLONE="pipelines/changeo-clone.sh"
ARG SHAZAM_THRESHOLD="pipelines/shazam-threshold.R"
ARG TIGGER_GENOTYPE="pipelines/tigger-genotype.R"
ARG PREPROCESS_PHIX="pipelines/preprocess-phix.sh"

# Install procotol data, utility scripts and pipelines
RUN PACKAGE="immcantation" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone -b default https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
	&& mv /tmp/${PACKAGE}/scripts/* /usr/local/bin/ \
    && mv /tmp/${PACKAGE}/${PRESTO_ABSEQ} /usr/local/bin/presto-abseq \
    && mv /tmp/${PACKAGE}/${CHANGEO_IGBLAST} /usr/local/bin/changeo-igblast \
    && mv /tmp/${PACKAGE}/${CHANGEO_CLONE} /usr/local/bin/changeo-clone \
    && mv /tmp/${PACKAGE}/${SHAZAM_THRESHOLD} /usr/local/bin/shazam-threshold \
    && mv /tmp/${PACKAGE}/${TIGGER_GENOTYPE} /usr/local/bin/tigger-genotype \
    && mv /tmp/${PACKAGE}/${PREPROCESS_PHIX} /usr/local/bin/preprocess-phix \
    && mv /tmp/${PACKAGE}/protocols /usr/local/share/protocols \
    && builds write -n ${PACKAGE} -v $(hg id -i /tmp/${PACKAGE}) \
    && rm -r /tmp/${PACKAGE}

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

# Install CD-HIT
RUN CDHIT=$(versions get -n cd-hit) \
    && wget -q --show-progress --no-check-certificate \
       https://github.com/weizhongli/cdhit/archive/V${CDHIT}.tar.gz \
    && tar -zxf V${CDHIT}.tar.gz \
    && (cd cdhit-${CDHIT} && make && make install) \
    && rm -r V${CDHIT}.tar.gz cdhit-${CDHIT}

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

# Install PHYLIP
RUN PHYLIP=$(versions get -n phylip) \
    && wget -q --show-progress --no-check-certificate \
       http://evolution.gs.washington.edu/phylip/download/phylip-${PHYLIP}.tar.gz \
    && tar -zxf phylip-${PHYLIP}.tar.gz \
    && (cd phylip-${PHYLIP}/src && make -f Makefile.unx install) \
    && mv phylip-${PHYLIP}/exe/* /usr/local/bin \
    && rm -r phylip-${PHYLIP}.tar.gz phylip-${PHYLIP}

# Install tbl2asn
RUN wget -q --show-progress --no-check-certificate \
    ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz \
    && gunzip linux64.tbl2asn.gz \
    && mv linux64.tbl2asn /usr/local/bin/tbl2asn \
    && chmod +x /usr/local/bin/tbl2asn

# Install AIRR reference libraries
RUN AIRR_PY=$(versions get -n airr-py) \
    && AIRR_R=$(versions get -n airr-r) \
    && pip3 install airr==${AIRR_PY} \
    && Rscript -e "devtools::install_version('airr', '${AIRR_R}', quiet=FALSE, clean=TRUE)"

# Install IgPhyML
RUN IGPHYML=$(versions get -n igphyml) \
    && wget -q --show-progress --no-check-certificate \
       https://bitbucket.org/kbhoehn/igphyml/get/${IGPHYML}.tar.gz -O igphyml-${IGPHYML}.tar.gz \
    && mkdir -p /usr/local/share/igphyml \
    && tar -zxf igphyml-${IGPHYML}.tar.gz -C /usr/local/share/igphyml --strip-components 1 \
    && (cd /usr/local/share/igphyml && ./make_phyml_blas_omp) \
    && rm -r igphyml-${IGPHYML}.tar.gz

# Install presto
RUN PACKAGE="presto" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone -b default https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && pip3 install --no-cache-dir /tmp/${PACKAGE} \
    && builds write -n ${PACKAGE} -v $(hg id -i /tmp/${PACKAGE}) \
    && rm -r /tmp/${PACKAGE}

# Install changeo
RUN PACKAGE="changeo" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone -b default https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && pip3 install --no-cache-dir /tmp/${PACKAGE} \
    && builds write -n ${PACKAGE} -v $(hg id -i /tmp/${PACKAGE}) \
    && rm -r /tmp/${PACKAGE}

# Install alakazam
RUN PACKAGE="alakazam" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone -b default https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && rinstall -p /tmp/${PACKAGE} \
    && builds write -n ${PACKAGE} -v $(hg id -i /tmp/${PACKAGE})

# Install shazam
RUN PACKAGE="shazam" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone -b default https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && rinstall -p /tmp/${PACKAGE} \
    && builds write -n ${PACKAGE} -v $(hg id -i /tmp/${PACKAGE})

# Install tigger
RUN PACKAGE="tigger" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone -b default https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && rinstall -p /tmp/${PACKAGE} \
    && builds write -n ${PACKAGE} -v $(hg id -i /tmp/${PACKAGE}) \
    && rm -r /tmp/${PACKAGE}

# Install rdi
RUN PACKAGE="rdi" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone https://bitbucket.org/cbolen1/rdicore /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && rinstall -p /tmp/${PACKAGE} \
    && builds write -n ${PACKAGE} -v $(hg id -i /tmp/${PACKAGE}) \
    && rm -r /tmp/${PACKAGE}

# Install RAbHIT
RUN PACKAGE="rabhit" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone --single-branch --branch master https://bitbucket.org/yaarilab/rabhit /tmp/${PACKAGE} \
    && rinstall -p /tmp/${PACKAGE} \
    && HERE=$(pwd) \
    && cd /tmp/${PACKAGE} \
    && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+) \
    && cd ${HERE} \
    && rm -r /tmp/${PACKAGE} 
    
# Install SCOPe
RUN PACKAGE="scoper" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone -b default https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && rinstall -p /tmp/${PACKAGE} \
    && builds write -n ${PACKAGE} -v $(hg id -i /tmp/${PACKAGE}) \
    && rm -r /tmp/${PACKAGE}

# Install prestor
RUN PACKAGE="prestor" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone -b default https://bitbucket.org/javh/prototype-prestor /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} \
    && rinstall -p /tmp/${PACKAGE} \
    && builds write -n ${PACKAGE} -v $(hg id -i /tmp/${PACKAGE}) \
    && rm -r /tmp/${PACKAGE}

# Download and build references databases
RUN fetch_igblastdb.sh -o /usr/local/share/igblast \
    && fetch_imgtdb.sh -o /usr/local/share/germlines/imgt \
    && imgt2igblast.sh -i /usr/local/share/germlines/imgt -o /usr/local/share/igblast \
    && fetch_phix.sh -o /usr/local/share/phix

# Setup environment
ENV PATH="${PATH}:/usr/local/share/igphyml/src"

# Setup users and permissions
RUN useradd magus -u 1000 -g users \
    && echo "magus ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/magus \
    && chmod 0440 /etc/sudoers.d/magus \
    && echo "umask 000" >> /home/magus/.bashrc
USER magus

# Set commands
CMD echo -e "Report version information:\n" \
            "  versions report\n" \
            "Report build details:\n" \
            "  builds report\n" \
            "Available pipeline commands:\n" \
            "  presto-abseq\n" \
            "  changeo-igblast\n" \
            "  changeo-clone\n" \
            "  tigger-genotype\n" \
            "  shazam-threshold\n" \
            "  preprocess-phix"
