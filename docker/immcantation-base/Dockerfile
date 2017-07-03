FROM fedora:25
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Dependencies, applications and scripts for the Immcantation framework."

# Versioning
COPY Version.yaml /Version.yaml
COPY tools/versions.py /usr/local/bin/versions

# Update and install required packages
RUN dnf -y update && dnf install -y \
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
	tar \
    texlive-scheme-basic \
	texlive-biblatex \
	texlive-titling \
	wget \
	&& dnf clean all

# Install utility scripts
RUN PACKAGE="immcantation" \
    && hg clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
	&& mv /tmp/${PACKAGE}/scripts/fastq2fasta.py /usr/local/bin/fastq2fasta.py \
	&& mv /tmp/${PACKAGE}/scripts/fetch_igblastdb.sh /usr/local/bin/fetch_igblastdb.sh \
	&& mv /tmp/${PACKAGE}/scripts/fetch_imgtdb.sh /usr/local/bin/fetch_imgtdb.sh \
	&& mv /tmp/${PACKAGE}/scripts/clean_imgtdb.py /usr/local/bin/clean_imgtdb.py \
	&& mv /tmp/${PACKAGE}/scripts/imgt2igblast.sh /usr/local/bin/imgt2igblast.sh \
	&& mv /tmp/${PACKAGE}/scripts/run_igblast.sh /usr/local/bin/run_igblast.sh \
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

# Download and build references databases
RUN fetch_igblastdb.sh -o /usr/local/share/igblast \
    && fetch_imgtdb.sh -o /usr/local/share/germlines/imgt \
    && imgt2igblast.sh -i /usr/local/share/germlines/imgt -o /usr/local/share/igblast
