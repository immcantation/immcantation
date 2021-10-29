FROM immcantation/suite:devel
LABEL maintainer="Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Immcantation Lab"
      
USER root        
# Remove this allele, for tigger to be able to find it
RUN sudo sed -i -e '/IGHV3-20\*04/{N;d;}' /usr/local/share/igblast/fasta/imgt_human_ig_v.fasta && \
    sed -i -e '/IGHV3-20\*04/{N;d;}' /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta && \
    imgt2igblast.sh -i /usr/local/share/germlines/imgt -o /usr/local/share/igblast
    
ENV NB_USER magus
ENV NB_UID 1000
ENV HOME /home/${NB_USER}
RUN mkdir ${HOME}/notebooks && \
    mkdir ${HOME}/data && \
    mkdir ${HOME}/results
VOLUME ${HOME}/results    
# Get notebook from repo
RUN PACKAGE="immcantation" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && mv /tmp/${PACKAGE}/training/*.ipynb ${HOME}/notebooks/. \
    && mv /tmp/${PACKAGE}/training/assets ${HOME}/notebooks/ \
    && rm -r /tmp/${PACKAGE}
        
RUN pip3 install --no-cache-dir rpy2 tzlocal notebook==6.* RISE && \
    pip3 install --no-cache-dir jupyter-contrib-nbextensions jupyter-nbextensions-configurator ipykernel bash_kernel && \
    python3 -m bash_kernel.install && \
    jupyter contrib nbextension install --sys-prefix && \
    jupyter nbextension install rise --py --sys-prefix && \
    jupyter nbextension enable splitcell/splitcell && \
    jupyter nbextensions_configurator enable

RUN Rscript -e "dir.create(path = Sys.getenv('R_LIBS_USER'), showWarnings = FALSE, recursive = TRUE)"
ARG R_DEPS="c('uuid', \
            'Seurat')"
RUN Rscript -e "install.packages(${R_DEPS}, lib=Sys.getenv('R_LIBS_USER'),clean=TRUE);devtools::install_github('RBigData/pbdZMQ');devtools::install_github(paste0('IRkernel/',c('repr', 'IRdisplay', 'IRkernel')));IRkernel::installspec(user = FALSE)"

WORKDIR ${HOME}/notebooks
EXPOSE 8888

COPY start-notebook.sh /usr/local/bin/
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
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
# Download dowsert tutorial data
RUN wget -v -O /home/${NB_USER}/immcantation-bcr_phylo_tutorial.zip -L https://yale.box.com/shared/static/96y1rlu4ujerprxr53qxw87g7ybc1y99 && \
    cd /home/${NB_USER}/ && \
    unzip -o immcantation-bcr_phylo_tutorial.zip && \
    rm -rf immcantation-bcr_phylo_tutorial.zip && \
    mv bcr_phylo_tutorial/* /home/${NB_USER}/data/. && \
    rm -rf bcr_phylo_tutorial  
#RUN cd ${HOME}/notebooks && \
#    jupyter nbconvert --ExecutePreprocessor.timeout=-1 --to notebook --execute intro-lab.ipynb && rm -rf results/*
CMD ["start-notebook.sh"]
