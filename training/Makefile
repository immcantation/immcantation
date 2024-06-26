# Convenience commands to run Immcanatation Lab + RStudio to
# work with the training notebooks
# Use from path: immcatation/training

IMMCANTATION_HOME=`realpath ..`
IMMCANTATION_IMAGE=immcantation/lab:devel
PASSWORD=immcantation

#------------
# Edit the local copies of the Rmd files
edit-notebooks:
	docker run --network host \
	-u $(id -u):$(id -g) \
	--name edit-notebooks \
	--detach \
	-v $(IMMCANTATION_HOME)/training:/home/magus/notebooks:z \
	--workdir /home/magus/notebooks \
	-e PASSWORD=$(PASSWORD) \
	--rm -p 8787:8787 $(IMMCANTATION_IMAGE)
# stop the edit-notebooks container
stop-edit-notebooks:
	docker ps | grep edit-notebooks | awk '{print $$NF}'| xargs docker stop
#-------------

#-------------
# Run rstudio
immcantation-rstudio:
	docker run --network host \
	-u $(id -u):$(id -g) \
	--name immcantation-rstudio \
	--detach \
	--workdir /home/magus/notebooks \
	-e PASSWORD=$(PASSWORD) \
	--rm -p 8787:8787 $(IMMCANTATION_IMAGE)

# stop the immcantation-rstudio container
stop-immcantation-rstudio:
	docker ps | grep immcantation-rstudio | awk '{print $$NF}'| xargs docker stop
#-------------

#-------------
# Render all local Rmd or rmd tutorials using 
# output_format='md_document'
render-tutorials:
	sudo rm -rf $(IMMCANTATION_HOME)/training/*_files/* && \
	docker run --network host \
	-u $(id -u):$(id -g) \
	--name render-tutorials \
	-v $(IMMCANTATION_HOME)/training:/home/magus/notebooks:z \
	--workdir /home/magus/notebooks \
	--rm $(IMMCANTATION_IMAGE) \
	Rscript -e "nbs <- list.files('.',pattern='[Rr]md$$',full.names=T);sapply(nbs, rmarkdown::render, output_format='md_document', params=list('eval_optional'=FALSE, 'cache'=FALSE), env=new.env())"

# stop the immcantation-rstudio container
stop-render-tutorials:
	docker ps | grep render-tutorials | awk '{print $$NF}'| xargs docker stop
#-------------
