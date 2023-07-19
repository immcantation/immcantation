ARG BASE_TAG=devel
ARG RS_USER=magus
ARG UID=1000
ARG GID=100
FROM immcantation/suite:${BASE_TAG}
LABEL maintainer="Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Immcantation + RStudio"

ARG BASE_TAG
ARG RS_USER
ARG UID
ARG GID
ENV RS_USER=$RS_USER

USER magus

RUN Rscript -e "dir.create(path = Sys.getenv('R_LIBS_USER'), showWarnings = FALSE, recursive = TRUE)"
RUN Rscript -e "devtools::install_github('rstudio/bookdown')"

###### RStudio #######
USER root

# Needed for RStudio server to open DESCRIPTION AND NAMESPACE
# Missing will yield error "File is binary rather than text so 
# cannot be opened by the source editor"
RUN echo "excludepkgs=packagename" /etc/dnf/dnf.conf
RUN dnf install -y \
            file \
            procps \
            pwgen \
            rstudio-server
            
ARG R_REPO="http://lib.stat.cmu.edu/R/CRAN/"
ARG R_DEPS="c('readxl', \
              'shiny')"
RUN echo "options(repos='${R_REPO}')" >> /usr/lib64/R/library/base/R/Rprofile \
    && mkdir -p /usr/share/doc/R/html \
    && Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"
  
## Make RStudio compatible with case when R is built from source
## (and thus is at /usr/local/bin/R), because RStudio doesn't obey
## path if a user apt-get installs a package    
RUN R_BIN=$(which R) \
    && mkdir -p /etc/rstudio \
    && echo "rsession-which-r=${R_BIN}" > /etc/rstudio/rserver.conf

#RUN cp /usr/lib/systemd/system/rstudio-server.service /etc/systemd/system/rstudio-server.service \
#   && sed -i '/Type=forking/a StartLimitIntervalSec=0' /etc/systemd/system/rstudio-server.service
# Create rstudio user
# useradd ${RS_USER} -u ${UID} --non-unique \
RUN PASSWD=$(pwgen 16 1) \
    && echo "${RS_USER}:${PASSWD}" | chpasswd \
    && usermod -G ${GID} ${RS_USER} \
    && echo "setwd('~/notebooks')" > /home/${RS_USER}/.Rprofile \
    && touch /etc/rstudio/rsession.conf \
    && echo "session-default-working-dir=~/notebooks" >> /etc/rstudio/rsession.conf \
    && echo "session-default-new-project-dir=~/notebooks" >>  /etc/rstudio/rsession.conf

# Fix for ERROR system error 13 (Permission denied) 
RUN mkdir -p /home/${RS_USER}/.local/share/rstudio \
    && chown -R ${RS_USER}:${GID} /home/${RS_USER}

COPY rstudio-server.sh /usr/bin/rstudio-server.sh
COPY run-rstudio.sh /usr/bin/run-rstudio.sh

RUN chmod +x /usr/bin/run-rstudio.sh
RUN chmod +x /usr/bin/rstudio-server.sh
RUN mkdir -p /var/log/rstudio-server
RUN export PATH=$PATH:/usr/local/bin

###### End RStudio #######

###### Tutorials #######

# Create databases without this allele, for tigger to be able to find it in the tutorial
RUN cp -r /usr/local/share/germlines/imgt /usr/local/share/germlines/imgt_test_tigger && \
    cp -r /usr/local/share/igblast /usr/local/share/igblast_test_tigger && \
    sudo sed -i -e '/IGHV3-20\*04/{N;d;}' /usr/local/share/igblast_test_tigger/fasta/imgt_human_ig_v.fasta && \
    sed -i -e '/IGHV3-20\*04/{N;d;}' /usr/local/share/germlines/imgt_test_tigger/human/vdj/imgt_human_IGHV.fasta && \
    imgt2igblast.sh -i /usr/local/share/germlines/imgt_test_tigger -o /usr/local/share/igblast_test_tigger

ENV HOME /home/${RS_USER}
ENV NB_USER ${RS_USER}
 
RUN mkdir ${HOME}/notebooks && \
    mkdir ${HOME}/data && \
    mkdir ${HOME}/results

# Get notebooks from repo
RUN PACKAGE="immcantation" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && mv /tmp/${PACKAGE}/training/*.Rmd ${HOME}/notebooks/. \
    && mv /tmp/${PACKAGE}/training/assets ${HOME}/notebooks/ \
    && rm -rf /tmp/${PACKAGE}

# Download intro lab data
RUN wget -v -O /home/${NB_USER}/immcantation-lab-intro.zip -L https://yale.box.com/shared/static/4bo611b70x8u92qvss1pypmcr9wmqil4 && \
    cd /home/${NB_USER}/ && \
    unzip -o immcantation-lab-intro.zip && \
    rm -rf immcantation-lab-intro.zip && \
    mv immcantation-lab-intro/example-input/* /home/${NB_USER}/data/. && \
    rm -rf immcantation-lab-intro
# Download BCR-Seurat tutorial data
RUN wget -v -O /home/${NB_USER}/immcantation-BCR-Seurat-tutorial.zip -L https://yale.box.com/shared/static/6s8gzxxpfb3vmwx2izn9qur9teultzvu.zip && \
    cd /home/${NB_USER}/ && \
    unzip -o immcantation-BCR-Seurat-tutorial.zip && \
    rm -rf immcantation-BCR-Seurat-tutorial.zip && \
    mv BCR.data.rds GEX.data.rds /home/${NB_USER}/data/. && \
    rm -rf immcantation-BCR-Seurat-tutorial
# Download 10x tutorial data
RUN wget -v -O /home/${NB_USER}/10x_data_2subj.zip -L http://clip.med.yale.edu/immcantation/examples/10x_data_2subj.zip && \
    cd /home/${NB_USER}/ && \
    unzip -o 10x_data_2subj.zip && \
    rm -rf 10x_data_2subj.zip && \
    mkdir -p  /home/${NB_USER}/data/10x_data_2subj && \
    mv sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_vdj_b_all_contig_igblast_db-pass.tsv \
    filtered_contig_annotations.csv \
    filtered_contig.fasta /home/${NB_USER}/data/10x_data_2subj/.
# Download dowsert tutorial data
RUN wget -v -O /home/${NB_USER}/immcantation-bcr_phylo_tutorial.zip -L https://yale.box.com/shared/static/96y1rlu4ujerprxr53qxw87g7ybc1y99 && \
    cd /home/${NB_USER}/ && \
    unzip -o immcantation-bcr_phylo_tutorial.zip && \
    rm -rf immcantation-bcr_phylo_tutorial.zip && \
    mv bcr_phylo_tutorial /home/${NB_USER}/data/. && \
    rm -rf __MACOSX

RUN chown -R ${RS_USER}:${GID} ${HOME}/notebooks \
    && chown -R ${RS_USER}:${GID} ${HOME}/data \
    && chown -R ${RS_USER}:${GID} ${HOME}/results

##### End Tutorials ########
VOLUME ${HOME}/results
WORKDIR ${HOME}/notebooks
ENV PASSWORD rstudio
CMD ["/usr/bin/run-rstudio.sh"]
