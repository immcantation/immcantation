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
#     && wget2 -q --force-progress --no-check-certificate \
#        http://evolution.gs.washington.edu/phylip/download/phylip-${PHYLIP}.tar.gz \
#     && tar -zxf phylip-${PHYLIP}.tar.gz \
#     && (cd phylip-${PHYLIP}/src && make -f Makefile.unx install CFLAGS='-fcommon') \
#     && mv phylip-${PHYLIP}/exe/* /usr/local/bin \
#     && rm -r phylip-${PHYLIP}.tar.gz phylip-${PHYLIP}

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
    ftp://ftp.ncbi.nlm.nih.gov/toolbox/ncbi_tools/converters/versions/2023-07-13/by_program/tbl2asn/linux64.tbl2asn.gz \
    && gunzip linux64.tbl2asn.gz \
    && mv linux64.tbl2asn /usr/local/bin/tbl2asn \
    && chmod +x /usr/local/bin/tbl2asn

# Install AIRR reference libraries
RUN AIRR_PY=$(versions get -n airr-py) \
    && AIRR_R=$(versions get -n airr-r) \
    && pip3 install airr==${AIRR_PY} \
    && Rscript -e "devtools::install_version('airr', '${AIRR_R}', quiet=FALSE, clean=TRUE)"

# Install alakazam
RUN ALAKAZAM_VER=$(versions get -n alakazam) \
    && Rscript -e "devtools::install_version('alakazam', '${ALAKAZAM_VER}', quiet=FALSE, clean=TRUE)"

# Install shazam
RUN SHAZAM_VER=$(versions get -n shazam) \
    && Rscript -e "devtools::install_version('shazam', '${SHAZAM_VER}', quiet=FALSE, clean=TRUE)"

# Install SCOPer
RUN SCOPER_VER=$(versions get -n scoper) \
    && Rscript -e "devtools::install_version('scoper', '${SCOPER_VER}', quiet=FALSE, clean=TRUE)"

# Install dowser
RUN DOWSER_VER=$(versions get -n dowser) \
    && Rscript -e "devtools::install_version('dowser', '${DOWSER_VER}', quiet=FALSE, clean=TRUE)"


# Install enchantr
RUN ENCHANTR_VER=$(versions get -n enchantr) \
    && Rscript -e "devtools::install_bitbucket('kleinstein/enchantr', '${ENCHANTR_VER}', quiet=FALSE, clean=TRUE)"

# Install IgPhyML
RUN IGPHYML=$(versions get -n igphyml) \
    && wget2 -q --force-progress --no-check-certificate \
       https://bitbucket.org/kleinstein/igphyml/get/${IGPHYML}.tar.gz -O igphyml-${IGPHYML}.tar.gz \
    && mkdir -p /usr/local/share/igphyml \
    && tar -zxf igphyml-${IGPHYML}.tar.gz -C /usr/local/share/igphyml --strip-components 1 \
    && (cd /usr/local/share/igphyml && ./make_phyml_blas_omp) \
    && rm -r igphyml-${IGPHYML}.tar.gz

# Setup environment
ENV PATH="${PATH}:/usr/local/share/igphyml/src"

# Set commands
CMD echo -e " Version information:\n"\
            "  versions report\n"\
            "Build stamps:\n"\
            "  builds report\n"\
            "Description of available pipelines:\n"\
            "  pipelines report"
