FROM kleinstein/immcantation:devel
LABEL maintainer="Susanna Marquez [susanna.marquez@yale.edu]" \
      description="Immcantation GUI"

USER root
      
ARG SHINY_URL="https://download3.rstudio.org/centos6.3/x86_64/shiny-server-1.5.7.907-rh6-x86_64.rpm"
ARG DOWNLOADS_DIR="/home/magus"


# Install R packages 
ARG R_DEPS="c('Rcpp', \
              'httpuv', \
              'shiny')"
ARG R_BUILD="library(devtools); install_deps(upgrade=FALSE, clean=TRUE); document(); install(build_vignettes=TRUE, clean=TRUE)"
            
RUN Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"

# Install shiny server
RUN mkdir -p /var/lib/shiny-server/bookmarks/shiny
RUN wget -O ${DOWNLOADS_DIR}/shiny-server.rpm ${SHINY_URL} \
    && dnf clean all \
    && dnf install -y --nogpgcheck ${DOWNLOADS_DIR}/shiny-server.rpm \
    && rm ${DOWNLOADS_DIR}/shiny-server.rpm

# Install R package immcatation-gui
RUN PACKAGE="immcantation-gui" \
    && rm -rf /tmp/${PACKAGE} \
    && hg clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && Rscript -e "setwd('/tmp/${PACKAGE}'); ${R_BUILD}" 
        
EXPOSE 3838

COPY shiny-server-app.R  /srv/shiny-server/immcantation/app.R
COPY shiny-server.sh  /usr/bin/shiny-server.sh
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

USER shiny

CMD ["/usr/bin/shiny-server.sh"]
