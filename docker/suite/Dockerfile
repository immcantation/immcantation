ARG BASE_TAG=devel
FROM immcantation/base:${BASE_TAG}
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Full suite of the Immcantation framework tools and dependencies."
ARG BASE_TAG
# Version
COPY Version.yaml /Version.yaml
COPY Pipeline.yaml /Pipeline.yaml
RUN builds write -n date -v "$(date +'%Y-%m-%d %T %Z')"

# Install procotol data, utility scripts and pipelines
RUN PACKAGE="immcantation" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && pipelines copy -s /tmp/${PACKAGE}/pipelines/ -t /usr/local/bin/ \
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

# Install muscle5
#RUN MUSCLE='5.1.0' \
#    && MUSCLE_SHORT=$(sed 's/\.0$//g' <<<"$MUSCLE") \
#    && wget2 -q --force-progress --no-check-certificate \
#       https://github.com/rcedgar/muscle/releases/download/${MUSCLE}/muscle${MUSCLE_SHORT}.linux_intel64 \
#    && chmod +x muscle${MUSCLE_SHORT}.linux_intel64 \
#	 && mv muscle${MUSCLE_SHORT}.linux_intel64 /usr/local/bin/muscle5

# Install vsearch
RUN VSEARCH=$(versions get -n vsearch) \
    && wget2 -q --force-progress --no-check-certificate \
       http://github.com/torognes/vsearch/releases/download/v${VSEARCH}/vsearch-${VSEARCH}-linux-x86_64.tar.gz \
    && tar -zxf vsearch-${VSEARCH}-linux-x86_64.tar.gz \
    && mv vsearch-${VSEARCH}-linux-x86_64/bin/vsearch /usr/local/bin/vsearch \
    && ln -s /usr/local/bin/vsearch /usr/local/bin/usearch \
    && rm -r vsearch-${VSEARCH}-linux-x86_64.tar.gz vsearch-${VSEARCH}-linux-x86_64

# Install CD-HIT
RUN CDHIT=$(versions get -n cd-hit) \
    && wget2 -q --force-progress --no-check-certificate \
       https://github.com/weizhongli/cdhit/archive/V${CDHIT}.tar.gz \
    && tar -zxf V${CDHIT}.tar.gz \
    && (cd cdhit-${CDHIT} && make && make install) \
    && rm -r V${CDHIT}.tar.gz cdhit-${CDHIT}

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

# # Install IgBLAST data
# RUN IGBLAST="1.19.0" \
#     && wget2 -q --force-progress --no-check-certificate \
#        ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST}/ncbi-igblast-${IGBLAST}-x64-linux.tar.gz \
#     && tar -zxf ncbi-igblast-${IGBLAST}-x64-linux.tar.gz \
#     && mkdir -p /usr/local/share/igblast \
#     && fetch_igblastdb.sh -o /usr/local/share/igblast \
#     && mv ncbi-igblast-${IGBLAST}/internal_data ncbi-igblast-${IGBLAST}/optional_file /usr/local/share/igblast \
#     && rm -r ncbi-igblast-${IGBLAST}-x64-linux.tar.gz ncbi-igblast-${IGBLAST}
#
# # Install IgBLAST executables
# COPY ncbi-igblast-1.19.0-fedora36-x64.zip /ncbi-igblast-1.19.0-fedora36-x64.zip
# RUN dnf install -y redhat-lsb-core sqlite-devel cmake
# RUN unzip -q /ncbi-igblast-1.19.0-fedora36-x64.zip \
#     && (cd ncbi-igblast-1.19.0-src/c++; make install) \
#     && rm -r ncbi-igblast-1.19.0-src \
#     && rm /ncbi-igblast-1.19.0-fedora36-x64.zip

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

# Install RAxML-NG
RUN RAXMLNG=$(versions get -n raxml-ng) \
    && wget2 -q --force-progress --no-check-certificate \
       https://github.com/amkozlov/raxml-ng/releases/download/1.2.0/raxml-ng_v${RAXMLNG}_linux_x86_64.zip \
    && unzip -q raxml-ng_v${RAXMLNG}_linux_x86_64.zip \
    && mv raxml-ng /usr/local/bin/raxml-ng \
    && chmod +x /usr/local/bin/raxml-ng \
    && rm -r raxml-ng_v${RAXMLNG}_linux_x86_64.zip

# Install tbl2asn, version 2023-07-13
# Important Announcement:
# Source: ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/README
# tbl2asn is no longer available for download. We encourage you
# to check out our updated version table2asn
# at https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/
# See its documentation at
# https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn_GFF/DOCUMENTATION/table2asn_readme.txt
RUN wget2 -q --force-progress --no-check-certificate \
    https://ftp.ncbi.nlm.nih.gov/toolbox/ncbi_tools/converters/versions/2023-07-13/by_program/tbl2asn/linux64.tbl2asn.gz \
    && gunzip linux64.tbl2asn.gz \
    && mv linux64.tbl2asn /usr/local/bin/tbl2asn \
    && chmod +x /usr/local/bin/tbl2asn

# Install AIRR reference libraries
RUN AIRR_PY=$(versions get -n airr-py) \
    && AIRR_R=$(versions get -n airr-r) \
    && pip3 install airr==${AIRR_PY} \
    && Rscript -e "devtools::install_version('airr', '${AIRR_R}', quiet=FALSE, clean=TRUE)"

# Install presto
RUN PACKAGE="presto" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://github.com/immcantation/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && pip3 install --no-cache-dir /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install changeo
RUN PACKAGE="changeo" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://github.com/immcantation/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && pip3 install --no-cache-dir /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install alakazam
RUN PACKAGE="alakazam" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://github.com/immcantation/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && rinstall -p /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install shazam
RUN PACKAGE="shazam" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && rinstall -p /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install tigger
RUN PACKAGE="tigger" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://github.com/immcantation/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && rinstall -p /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install SCOPer
RUN PACKAGE="scoper" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && rinstall -p /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install dowser
RUN PACKAGE="dowser" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && rinstall -p /tmp/${PACKAGE} -u never \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install enchantr
RUN PACKAGE="enchantr" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://github.com/immcantation/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && rinstall -p /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install prestor
# Note: Package ‘captioner’ was removed from the CRAN repository.
# Archived on 2023-08-18
RUN PACKAGE="prestor" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://github.com/immcantation/prestor /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && Rscript -e "devtools::install_github('adletaw/captioner')" \
    && rinstall -p /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install IgPhyML
RUN IGPHYML=$(versions get -n igphyml) \
    && wget2 -q --force-progress --no-check-certificate \
       https://bitbucket.org/kleinstein/igphyml/get/${IGPHYML}.tar.gz -O igphyml-${IGPHYML}.tar.gz \
    && mkdir -p /usr/local/share/igphyml \
    && tar -zxf igphyml-${IGPHYML}.tar.gz -C /usr/local/share/igphyml --strip-components 1 \
    && (cd /usr/local/share/igphyml && ./make_phyml_blas_omp) \
    && rm -r igphyml-${IGPHYML}.tar.gz

# Install rdi
RUN PACKAGE="rdi" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://bitbucket.org/cbolen1/rdicore /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && rinstall -p /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install RAbHIT
RUN RABHIT=$(versions get -n rabhit) \
    && Rscript -e "devtools::install_version('rabhit', '${RABHIT}', quiet=FALSE, clean=TRUE)"

# Install piglet
RUN PIGLET=$(versions get -n piglet) \
    && Rscript -e "devtools::install_version('piglet', '${PIGLET}', quiet=FALSE, clean=TRUE)"

# Install Seurat
RUN SEURAT=$(versions get -n seurat) \
    && Rscript -e "devtools::install_version('Seurat', '${SEURAT}', quiet=FALSE, clean=TRUE, upgrade='never')"

# Download and build references databases
RUN fetch_imgtdb.sh -o /usr/local/share/germlines/imgt \
    && imgt2igblast.sh -i /usr/local/share/germlines/imgt -o /usr/local/share/igblast \
    && fetch_phix.sh -o /usr/local/share/phix

# Patch IgBLAST database
#RUN wget2 -q --force-progress --no-check-certificate \
#    ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/patch/optional_file/mouse_gl.aux \
#    && mv -f mouse_gl.aux /usr/local/share/igblast/optional_file

# Setup environment
ENV PATH="${PATH}:/usr/local/share/igphyml/src"

# Setup users and permissions
RUN useradd magus -u 1000 -g users \
    && echo "magus ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/magus \
    && chmod 0440 /etc/sudoers.d/magus \
    && echo "umask 000" >> /home/magus/.bashrc \
    && echo "umask 000" >> /root/.bashrc
# USER magus

# Set commands
CMD echo -e " Version information:\n"\
            "  versions report\n"\
            "Build stamps:\n"\
            "  builds report\n"\
            "Description of available pipelines:\n"\
            "  pipelines report"
