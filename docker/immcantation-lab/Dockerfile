FROM kleinstein/immcantation:devel
LABEL maintainer="Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Immcantation Lab"
      
USER root        

ENV NB_USER magus
ENV NB_UID 1000
ENV HOME /home/${NB_USER}
RUN mkdir ${HOME}/notebooks && \
    mkdir ${HOME}/data && \
    mkdir ${HOME}/results
VOLUME ${HOME}/results    
# Get notebook from repo
# TODO: clone from repo
#COPY *.ipynb ${HOME}/notebooks/
# Install procotol data, utility scripts and pipelines
RUN PACKAGE="immcantation" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone -b default https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && mv /tmp/${PACKAGE}/training/intro-lab.ipynb ${HOME}/notebooks \
    && mv /tmp/${PACKAGE}/training/assets ${HOME}/notebooks/ \
    && rm -r /tmp/${PACKAGE}
    
USER root
WORKDIR ${HOME}/notebooks
RUN pip3 install --no-cache-dir rpy2 tzlocal notebook==5.* RISE jupyter-contrib-nbextensions jupyter-nbextensions-configurator ipykernel bash_kernel && \
    python3 -m bash_kernel.install && \
    jupyter contrib nbextension install --sys-prefix && \
    jupyter nbextension install rise --py --sys-prefix && \
    jupyter nbextension enable splitcell/splitcell && \
    jupyter nbextensions_configurator enable
RUN jupyter notebook --generate-config -y
RUN Rscript -e "dir.create(path = Sys.getenv('R_LIBS_USER'), showWarnings = FALSE, recursive = TRUE)"
ARG R_DEPS="c('repr', \
            'IRdisplay', \
            'evaluate', \
            'crayon', \
            'pbdZMQ', \
            'devtools', \
            'uuid', \
            'digest')"
RUN Rscript -e "install.packages(${R_DEPS}, lib=Sys.getenv('R_LIBS_USER'),clean=TRUE);devtools::install_github('IRkernel/IRkernel');IRkernel::installspec(user = FALSE)"
EXPOSE 8888
COPY start-notebook.sh /usr/local/bin/
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
RUN wget -v -O /home/${NB_USER}/immcantation-webinar.zip -L https://yale.box.com/shared/static/4bo611b70x8u92qvss1pypmcr9wmqil4.zip && \
    cd /home/${NB_USER}/ && \
    unzip -o immcantation-webinar.zip && \
    rm -rf immcantation-webinar.zip && \
    mv immcantation-webinar/example-input/* /home/${NB_USER}/data/. && \
    rm -rf immcantation-webinar
#RUN cd ${HOME}/notebooks && \
#    jupyter nbconvert --ExecutePreprocessor.timeout=-1 --to notebook --execute intro-lab.ipynb && \
#    rm -rf results/*
CMD ["start-notebook.sh"]
