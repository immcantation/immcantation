FROM fedora:25
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Dependencies, applications and scripts for the Immcantation framework."

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
    findutils \
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
