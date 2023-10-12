ARG BASE_TAG=devel
FROM immcantation/base:${BASE_TAG}
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Immcantation framework tools needed for airrflow."
ARG BASE_TAG
# Version
COPY Version.yaml /Version.yaml
RUN builds write -n date -v "$(date +'%Y-%m-%d %T %Z')"

# # Install PHYLIP
# RUN PHYLIP=$(versions get -n phylip) \
#     && wget -q --show-progress --no-check-certificate \
#        http://evolution.gs.washington.edu/phylip/download/phylip-${PHYLIP}.tar.gz \
#     && tar -zxf phylip-${PHYLIP}.tar.gz \
#     && (cd phylip-${PHYLIP}/src && make -f Makefile.unx install CFLAGS='-fcommon') \
#     && mv phylip-${PHYLIP}/exe/* /usr/local/bin \
#     && rm -r phylip-${PHYLIP}.tar.gz phylip-${PHYLIP}

# Install RAxML-NG
RUN RAXMLNG=$(versions get -n raxml-ng) \
    && wget -q --show-progress --no-check-certificate \
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
RUN wget -q --show-progress --no-check-certificate \
    ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/versions/2023-07-13/by_program/tbl2asn/linux64.tbl2asn.gz \
    && gunzip linux64.tbl2asn.gz \
    && mv linux64.tbl2asn /usr/local/bin/tbl2asn \
    && chmod +x /usr/local/bin/tbl2asn

# Install AIRR reference libraries
RUN AIRR_PY=$(versions get -n airr-py) \
    && AIRR_R=$(versions get -n airr-r) \
    && pip3 install airr==${AIRR_PY} \
    && Rscript -e "devtools::install_version('airr', '${AIRR_R}', quiet=FALSE, clean=TRUE)"

# Install alakazam
RUN PACKAGE="alakazam" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
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
    && git clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && rinstall -p /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install prestor
# Note: Package ‘captioner’ was removed from the CRAN repository.
# Archived on 2023-08-18
RUN PACKAGE="prestor" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://bitbucket.org/kleinstein/prestor /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && Rscript -e "devtools::install_github('adletaw/captioner')" \
    && rinstall -p /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Install IgPhyML
RUN IGPHYML=$(versions get -n igphyml) \
    && wget -q --show-progress --no-check-certificate \
       https://bitbucket.org/kleinstein/igphyml/get/${IGPHYML}.tar.gz -O igphyml-${IGPHYML}.tar.gz \
    && mkdir -p /usr/local/share/igphyml \
    && tar -zxf igphyml-${IGPHYML}.tar.gz -C /usr/local/share/igphyml --strip-components 1 \
    && (cd /usr/local/share/igphyml && ./make_phyml_blas_omp) \
    && rm -r igphyml-${IGPHYML}.tar.gz

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
