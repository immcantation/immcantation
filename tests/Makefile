.PHONY: presto changeo pipelines nfpipelines nfpipelines-resume nfpipelines-clean stop

# Convenience commands to test Immcantation release

# Run  tests without warnings
# make pipelines PYTHON_WARNINGS=ignore
# make pipelines  PYTHON_WARNINGS="ignore:::airr"
# make pipelines BATS_TEST="changeo-igblast" PYTHON_WARNINGS="ignore"
# make pipelines BATS_TEST="changeo-10x" PYTHON_WARNINGS="ignore:::airr"
# make pipelines BATS_TEST="changeo-igblast|changeo-clone" PYTHON_WARNINGS="ignore"

IMMCANTATION_HOME=`realpath ..`
IMMCANTATION_IMAGE ?= immcantation/suite:devel
PYTHON_WARNINGS ?= default
BATS_TEST ?= .

#------------
all: clean presto changeo pipelines

# make igblast PYTHON_WARNINGS=ignore
igblast:
	make pipelines BATS_TEST="changeo-igblast" PYTHONWARNINGS="$(PYTHON_WARNINGS)" IMAGE="$(IMMCANTATION_IMAGE)"

# make clone PYTHON_WARNINGS=ignore
clone:
	make pipelines BATS_TEST="changeo-clone" PYTHONWARNINGS="$(PYTHON_WARNINGS)" IMAGE="$(IMMCANTATION_IMAGE)"

# make tigger PYTHON_WARNINGS=ignore
tigger:
	make pipelines BATS_TEST="tigger-genotype" PYTHONWARNINGS="$(PYTHON_WARNINGS)" IMAGE="$(IMMCANTATION_IMAGE)"


clean:
	rm -rf run/*
#------------
# test_presto.sh
#------------
presto:
	docker run --network host \
	-u $(id -u):$(id -g) \
	--name prerelease-tests \
	-v $(IMMCANTATION_HOME)/tests:/data:z \
	--workdir /data \
	-e PYTHONWARNINGS=$(PYTHON_WARNINGS) \
	--rm $(IMMCANTATION_IMAGE) /bin/bash -c "\
	dnf install -y bats \
	&& ./test_presto.sh"

presto-gzip:
	docker run --network host \
	-u $(id -u):$(id -g) \
	--name prerelease-tests \
	-v $(IMMCANTATION_HOME)/tests:/data:z \
	--workdir /data \
	-e PYTHONWARNINGS=$(PYTHON_WARNINGS) \
	--rm $(IMMCANTATION_IMAGE) /bin/bash -c "\
	dnf install -y bats \
	&& ./test_presto-gzip.sh"

#------------
# test_changeo-airr.sh
#------------
changeo:
	docker run --network host \
	-u $(id -u):$(id -g) \
	--name prerelease-tests \
	-v $(IMMCANTATION_HOME)/tests:/data:z \
	--workdir /data \
	-e PYTHONWARNINGS=$(PYTHON_WARNINGS) \
	--rm $(IMMCANTATION_IMAGE) /bin/bash -c "\
	dnf install -y bats && ./test_changeo-airr.sh"

#------------
# test docker container and pipelines
#------------
pipelines:
	PYTHONWARNINGS="$(PYTHON_WARNINGS)" IMAGE="$(IMMCANTATION_IMAGE)" bats --jobs 1 --filter "$(BATS_TEST)" test_pipelines.sh

# Nextflow pipeline tests
nfpipelines:
	DATA_DIR=$$(realpath data) RUN_DIR=$$(realpath run) \
	nextflow run test_nfpipelines.nf \
		-profile docker \
		--image "$(IMMCANTATION_IMAGE)" \
		--pythonwarnings "$(PYTHON_WARNINGS)" \
		-c <(echo "docker.runOptions = '-u \$$(id -u):\$$(id -g) -v \$${DATA_DIR}:/data:z -v \$${RUN_DIR}:/scratch:z'")

# Resume nextflow tests
nfpipelines-resume:
	DATA_DIR=$$(realpath data) RUN_DIR=$$(realpath run) \
	nextflow run test_nfpipelines.nf \
		-profile docker \
		-resume \
		--image "$(IMMCANTATION_IMAGE)" \
		--pythonwarnings "$(PYTHON_WARNINGS)" \
		-c <(echo "docker.runOptions = '-u \$$(id -u):\$$(id -g) -v \$${DATA_DIR}:/data:z -v \$${RUN_DIR}:/scratch:z'")

# Clean nextflow work directory
nfpipelines-clean:
	rm -rf work/ .nextflow* .nextflow.log*

# stop the prerelease-tests container
stop:
	docker ps | grep prerelease-tests | awk '{print $$NF}'| xargs docker stop
