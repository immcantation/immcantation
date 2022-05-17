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
    && rm -rf /tmp/${PACKAGE}

RUN pip3 install --no-cache-dir "nbconvert<6" "jupyter-client<7"
RUN pip3 install --no-cache-dir "jinja2<3.1.0"
RUN pip3 install --no-cache-dir rpy2 tzlocal notebook RISE
RUN pip3 install --no-cache-dir pip install jupyter_contrib_nbextensions \
    jupyter_nbextensions_configurator \
    ipykernel bash_kernel
RUN python3 -m bash_kernel.install
RUN jupyter contrib nbextension install --sys-prefix
# template_path workaround https://github.com/ipython-contrib/jupyter_contrib_nbextensions/issues/1529
# RUN PYTHON_V=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
RUN sed -i 's/template_path/template_paths/g' /usr/local/lib/python3.9/site-packages/jupyter_contrib_nbextensions/nbconvert_support/toc2.py && \
    sed -i 's/template_path/template_paths/g' /usr/local/lib/python3.9/site-packages/jupyter_contrib_nbextensions/nbconvert_support/exporter_inliner.py

RUN dnf install -y R-rgeos
RUN Rscript -e "dir.create(path = Sys.getenv('R_LIBS_USER'), showWarnings = FALSE, recursive = TRUE)"
ARG R_DEPS="c('uuid', \
            'Seurat')"
RUN Rscript -e "install.packages(${R_DEPS}, lib=Sys.getenv('R_LIBS_USER'),clean=TRUE);devtools::install_github('RBigData/pbdZMQ');devtools::install_github(paste0('IRkernel/',c('repr', 'IRdisplay', 'IRkernel')));IRkernel::installspec(user = FALSE)"

WORKDIR ${HOME}/notebooks
EXPOSE 8888

COPY start-notebook.sh /usr/local/bin/
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

RUN jupyter nbextension enable toc2/toc2
RUN jupyter nbextension enable splitcell/splitcell
# RUN jupyter nbextensions_configurator enable --sys-prefix

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
#RUN cd ${HOME}/notebooks && \
#    jupyter nbconvert --ExecutePreprocessor.timeout=-1 --to notebook --execute intro-lab.ipynb && rm -rf results/*
CMD ["start-notebook.sh"]
