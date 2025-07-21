
# Immcantation Docker Makefile

.PHONY: build clean release help
help:
	@echo "=================================="
	@echo "   Immcantation Docker Makefile   "
	@echo "=================================="
	@echo ""
	@echo "This Makefile is used to build Docker images for the Immcantation suite."
	@echo "It provides:"
	@echo "  - a generic build target for any supported image"
	@echo "  - a release target to build all images"
	@echo "  - a clean target to delete the log files"
	@echo ""
	@echo "Usage: make build IMAGE=<image> [TAG=<tag>] [NO_CACHE=1]"
	@echo "Default values: TAG=devel, NO_CACHE=0 (build uses cache by default)"
	@echo "Possible IMAGE values: base, suite, lab, test, airrflow"
	@echo ""
	@echo "Examples:"
	@echo "  make build IMAGE=base TAG=devel NO_CACHE=1"
	@echo "  make clean"

TAG ?= devel
NO_CACHE_ARG =
ifeq ($(NO_CACHE),1)
  NO_CACHE_ARG = --no-cache
endif

build:
	@if [ -z "$(IMAGE)" ]; then \
	  echo "Error: IMAGE variable not set. Usage: make build IMAGE=base [TAG=...] [NO_CACHE=1]"; \
	  exit 1; \
	fi
	@case " $(IMAGE) " in \
	  " base "|" suite "|" lab "|" test "|" airrflow ") ;; \
	  *) echo "Error: IMAGE must be one of: base, suite, lab, test, airrflow"; exit 1 ;; \
	esac
	@docker build $(IMAGE)/. -t immcantation/$(IMAGE):$(TAG) $(NO_CACHE_ARG) --progress=plain -f $(IMAGE)/Dockerfile 2>&1 | tee $(IMAGE)-$(TAG).log

release:
	$(foreach img,base suite lab test airrflow,$(MAKE) build IMAGE=$(img);)

clean:
	rm -f base-*.log suite-*.log lab-*.log test-*.log airrflow-*.log
