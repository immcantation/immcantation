.PHONY: presto changeo stop

# Convenience commands to test Immcantation release

IMMCANTATION_HOME=`realpath ..`
IMMCANTATION_IMAGE=immcantation/suite:4.6.0

#------------
all: presto changeo

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
	--rm $(IMMCANTATION_IMAGE) /bin/bash -c "\
	dnf install -y bats && ./test_presto.sh"

#------------
# test_changeo-airr.sh
#------------
changeo:
	docker run --network host \
	-u $(id -u):$(id -g) \
	--name prerelease-tests \
	-v $(IMMCANTATION_HOME)/tests:/data:z \
	--workdir /data \
	--rm $(IMMCANTATION_IMAGE) /bin/bash -c "\
	dnf install -y bats && ./test_changeo-airr.sh"


# stop the immcantation-rstudio container
stop:
	docker ps | grep prerelease-tests | awk '{print $$NF}'| xargs docker stop
