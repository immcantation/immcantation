.PHONY: help edit-notebooks stop-edit-notebooks immcantation-rstudio stop-immcantation-rstudio clean render-tutorials stop-render-tutorials check-docker

# Immcantation Lab + RStudio Training Environment
# Use from path: immcantation/training
#
# Common usage:
#   make help                    - Show this help
#   make edit-notebooks          - Start RStudio to edit notebooks
#   make immcantation-rstudio    - Start RStudio with full environment
#   make render-tutorials        - Render all tutorials to markdown
#   make clean                   - Clean generated files
#   make render-tutorials IMMCANTATION_IMAGE=immcantation/lab:4.7.0 FILES='intro-lab.Rmd'

# Configuration - Override these as needed
IMMCANTATION_HOME ?= $(realpath ..)
IMMCANTATION_IMAGE ?= immcantation/lab:devel
PASSWORD ?= immcantation
PORT ?= 8787
USER_ID ?= $(id -u)
GROUP_ID ?= $(id -g)

# Colors for output
GREEN := \033[0;32m
YELLOW := \033[0;33m
RED := \033[0;31m
NC := \033[0m # No Color

#------------
# Help target - default when running 'make'
help:
	@echo -e "$(GREEN)Immcantation Training Environment$(NC)"
	@echo ""
	@echo -e "$(YELLOW)Available targets:$(NC)"
	@echo "  help                    - Show this help message"
	@echo "  edit-notebooks          - Start RStudio to edit local notebook files"
	@echo "  immcantation-rstudio    - Start RStudio with full Immcantation environment"
	@echo "  render-tutorials        - Render all R Markdown tutorials to markdown"
	@echo "                           Use FILES='file1.Rmd file2.Rmd' to render specific files"
	@echo "  clean                   - Remove generated tutorial files and cache"
	@echo "  stop-edit-notebooks     - Stop the edit-notebooks container"
	@echo "  stop-immcantation-rstudio - Stop the immcantation-rstudio container"
	@echo "  stop-render-tutorials   - Stop the render-tutorials container"
	@echo ""
	@echo -e "$(YELLOW)Configuration (override with make VAR=value):$(NC)"
	@echo "  IMMCANTATION_HOME=$(IMMCANTATION_HOME)"
	@echo "  IMMCANTATION_IMAGE=$(IMMCANTATION_IMAGE)"
	@echo "  PASSWORD=$(PASSWORD)"
	@echo "  PORT=$(PORT)"
	@echo "  USER_ID=$(USER_ID)"
	@echo "  GROUP_ID=$(GROUP_ID)"
	@echo ""
	@echo -e "$(YELLOW)Common Examples:$(NC)"
	@echo "  make edit-notebooks                              # Start RStudio with default settings"
	@echo "  make edit-notebooks PASSWORD=mypassword         # Use custom password"
	@echo "  make immcantation-rstudio PORT=8888              # Use different port"
	@echo "  make render-tutorials                            # Render all tutorials"
	@echo "  make render-tutorials FILES='intro-lab.Rmd'     # Render one tutorial"
	@echo "  make render-tutorials FILES='intro-lab.Rmd 10x_tutorial.Rmd'  # Render multiple"
	@echo "  make clean                                       # Clean generated files"
	@echo ""
	@echo -e "$(YELLOW)Quick Start:$(NC)"
	@echo "  1. Run 'make edit-notebooks' to start RStudio"
	@echo "  2. Open browser to http://localhost:$(PORT)"
	@echo "  3. Login with username: magus, password: $(PASSWORD)"
	@echo "  4. Edit .Rmd files in the notebooks directory"
	@echo "  5. Run 'make render-tutorials' to generate markdown files"
	@echo ""
	@echo -e "$(YELLOW)Troubleshooting:$(NC)"
	@echo "  - If port $(PORT) is busy: make edit-notebooks PORT=8888"
	@echo "  - If permission errors: check USER_ID/GROUP_ID match your system"
	@echo "  - To stop containers: make stop-edit-notebooks or make stop-immcantation-rstudio"

#------------
# Check if Docker is available
check-docker:
	@which docker > /dev/null || (echo -e "$(RED)Error: Docker not found. Please install Docker first.$(NC)" && exit 1)
	@docker info > /dev/null 2>&1 || (echo -e "$(RED)Error: Docker daemon not running. Please start Docker.$(NC)" && exit 1)

#------------
# Delete rendered tutorial files and cache
clean:
	@echo -e "$(YELLOW)Cleaning generated files...$(NC)"
	@sudo rm -rf $(IMMCANTATION_HOME)/training/*_files $(IMMCANTATION_HOME)/training/*_cache
	@echo -e "$(GREEN)Clean complete.$(NC)"

#------------
# Edit the local copies of the Rmd files
edit-notebooks: check-docker
	@echo -e "$(YELLOW)Starting RStudio for notebook editing...$(NC)"
	@echo -e "$(GREEN)RStudio will be available at: http://localhost:$(PORT)$(NC)"
	@echo -e "$(GREEN)Username: magus, Password: $(PASSWORD)$(NC)"
	@docker run --network host \
		-u $(USER_ID):$(GROUP_ID) \
		--name edit-notebooks \
		--detach \
		-v $(IMMCANTATION_HOME)/training:/home/magus/notebooks:z \
		--workdir /home/magus/notebooks \
		-e PASSWORD=$(PASSWORD) \
		--rm -p $(PORT):8787 $(IMMCANTATION_IMAGE) || \
		(echo -e "$(RED)Error: Failed to start container. Check if port $(PORT) is already in use.$(NC)" && exit 1)
	@echo -e "$(GREEN)Container 'edit-notebooks' started successfully.$(NC)"

# Stop the edit-notebooks container
stop-edit-notebooks:
	@echo -e "$(YELLOW)Stopping edit-notebooks container...$(NC)"
	@docker ps --filter "name=edit-notebooks" --format "{{.Names}}" | xargs -r docker stop
	@echo -e "$(GREEN)Container stopped.$(NC)"

#-------------
# Run RStudio with full Immcantation environment
immcantation-rstudio: check-docker
	@echo -e "$(YELLOW)Starting Immcantation RStudio environment...$(NC)"
	@echo -e "$(GREEN)RStudio will be available at: http://localhost:$(PORT)$(NC)"
	@echo -e "$(GREEN)Username: magus, Password: $(PASSWORD)$(NC)"
	@docker run --network host \
		-u $(USER_ID):$(GROUP_ID) \
		--name immcantation-rstudio \
		--detach \
		--workdir /home/magus/notebooks \
		-e PASSWORD=$(PASSWORD) \
		--rm -p $(PORT):8787 $(IMMCANTATION_IMAGE) || \
		(echo -e "$(RED)Error: Failed to start container. Check if port $(PORT) is already in use.$(NC)" && exit 1)
	@echo -e "$(GREEN)Container 'immcantation-rstudio' started successfully.$(NC)"

# Stop the immcantation-rstudio container
stop-immcantation-rstudio:
	@echo -e "$(YELLOW)Stopping immcantation-rstudio container...$(NC)"
	@docker ps --filter "name=immcantation-rstudio" --format "{{.Names}}" | xargs -r docker stop
	@echo -e "$(GREEN)Container stopped.$(NC)"

#-------------
# Render R Markdown tutorials - all files or specific files if FILES is provided
# make render-tutorials FILES='tutorial1.Rmd tutorial2.Rmd tutorial3.Rmd'
render-tutorials: check-docker
ifdef FILES
	@echo -e "$(YELLOW)Rendering specific tutorials: $(FILES)$(NC)"
	@for file in $(FILES); do \
		if [ ! -f "$(IMMCANTATION_HOME)/training/$$file" ]; then \
			echo -e "$(RED)Error: File $$file not found in training directory$(NC)"; \
			exit 1; \
		fi; \
		base=$$(basename $$file .Rmd); base=$$(basename $$base .rmd); \
		echo -e "$(YELLOW)Cleaning output for $$file...$(NC)"; \
		sudo rm -rf $(IMMCANTATION_HOME)/training/$${base}_files $(IMMCANTATION_HOME)/training/$${base}_cache 2>/dev/null || true; \
	done
	@docker run --network host \
		-u $(USER_ID):$(GROUP_ID) \
		--name render-tutorials \
		-v $(IMMCANTATION_HOME)/training:/home/magus/notebooks:z \
		--workdir /home/magus/notebooks \
		--rm $(IMMCANTATION_IMAGE) \
		Rscript -e "files <- c($(shell echo '$(FILES)' | sed "s/[^ ]*/'\0'/g" | tr ' ' ',')); for(f in files) { if(file.exists(f)) { cat('Rendering:', f, '\n'); rmarkdown::render(f, output_format='md_document', params=list('eval_optional'=FALSE, 'cache'=FALSE), env=new.env()) } else { cat('File not found:', f, '\n') } }" || \
		(echo -e "$(RED)Error: Rendering failed.$(NC)" && exit 1)
	@echo -e "$(GREEN)Specified tutorials rendered successfully.$(NC)"
else
	@echo -e "$(YELLOW)Rendering all tutorials to markdown...$(NC)"
	@sudo rm -rf $(IMMCANTATION_HOME)/training/*_files $(IMMCANTATION_HOME)/training/*_cache 2>/dev/null || true
	@docker run --network host \
		-u $(USER_ID):$(GROUP_ID) \
		--name render-tutorials \
		-v $(IMMCANTATION_HOME)/training:/home/magus/notebooks:z \
		--workdir /home/magus/notebooks \
		--rm $(IMMCANTATION_IMAGE) \
		Rscript -e "nbs <- list.files('.',pattern='[Rr]md$$',full.names=T); if(length(nbs)==0) {cat('No R Markdown files found\n'); quit()} else {cat('Found', length(nbs), 'files to render\n'); sapply(nbs, function(x) {cat('Rendering:', x, '\n'); rmarkdown::render(x, output_format='md_document', params=list('eval_optional'=FALSE, 'cache'=FALSE), env=new.env())})}" || \
		(echo -e "$(RED)Error: Rendering failed.$(NC)" && exit 1)
	@echo -e "$(GREEN)All tutorials rendered successfully.$(NC)"
endif

# Default target
.DEFAULT_GOAL := help