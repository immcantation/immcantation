IMMCANTATION_HOME=`realpath ..`
IMMCANTATION_IMAGE=immcantation/lab:devel

edit-notebooks:
	docker run --network=host \
	--name edit-notebooks \
	--detach \
	-v $(IMMCANTATION_HOME)/training:/home/magus/notebooks:z \
	--workdir /home/magus/notebooks \
	-e PASSWORD=immcantation \
	--rm -p 8787:8787 $(IMMCANTATION_IMAGE)

stop-edit-notebooks:
	docker ps | grep edit-notebooks | awk '{print $$NF}'| xargs docker stop

render-tutorials:
	docker run \
	--name tutorials \
	-v $(IMMCANTATION_HOME)/training:/home/magus/notebooks:z \
	--workdir /home/magus/notebooks \
	--rm $(IMMCANTATION_IMAGE) \
	Rscript -e "nbs <- list.files('.',pattern='Rmd$$',full.names=T);sapply(nbs, rmarkdown::render, output_format='md_document', env=new.env())"


render-10x_tutorial:
	sudo rm -rf $(IMMCANTATION_HOME)/training/10x_tutorial_files/* && \
	docker run \
	--name 10x_tutorial \
	-v $(IMMCANTATION_HOME)/training:/home/magus/notebooks:z \
	--workdir /home/magus/notebooks \
	--rm $(IMMCANTATION_IMAGE) \
	Rscript -e "nbs <- '10x_tutorial.Rmd';sapply(nbs, rmarkdown::render, output_format='md_document', params=list('eval_optional'=FALSE, 'cache'=FALSE), env=new.env())"

render-BCR_Seurat_tutorial:
	docker run \
	--name 10x_tutorial \
	-v $(IMMCANTATION_HOME)/training:/home/magus/notebooks:z \
	--workdir /home/magus/notebooks \
	--rm $(IMMCANTATION_IMAGE) \
	Rscript -e "nbs <- 'BCR_Seurat_tutorial.Rmd';sapply(nbs, rmarkdown::render, output_format='md_document', params=list('eval_optional'=FALSE, 'cache'=FALSE), env=new.env())"

render-intro_tutorial:
	sudo rm -rf $(IMMCANTATION_HOME)/training/intro-lab_files/* && \
	docker run \
	--name intro_tutorial \
	-v $(IMMCANTATION_HOME)/training:/home/magus/notebooks:z  \
	--workdir /home/magus/notebooks \
	--rm $(IMMCANTATION_IMAGE) \
	Rscript -e "nbs <- 'intro-lab.Rmd';sapply(nbs, rmarkdown::render, output_format='md_document', params=list('eval_optional'=TRUE, 'cache'=FALSE), env=new.env())"

