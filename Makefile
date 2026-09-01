# PureCLIP — developer convenience targets.
# The build itself lives in src/CMakeLists.txt; this is a thin wrapper so that
# `make test` works from a clean checkout. See TESTING.md.

BUILD   ?= build
JOBS    ?= 8
CMAKE   ?= cmake
CTEST   ?= ctest

.PHONY: all build test test-all test-truth golden clean

all: build

build:
	@$(CMAKE) -S src -B $(BUILD) -DCMAKE_BUILD_TYPE=Release
	@$(CMAKE) --build $(BUILD) -j $(JOBS)

## Fast correctness check on the synthetic fixture (~seconds).
test: build
	@$(CTEST) --test-dir $(BUILD) -L tier1 --output-on-failure

## Everything, including the real chrM data (~30 s).
test-all: build
	@$(CTEST) --test-dir $(BUILD) --output-on-failure

## Only the ground-truth check: are the implanted crosslinks recovered?
test-truth: build
	@$(CTEST) --test-dir $(BUILD) -L truth --output-on-failure

## Re-record reference output. Only after confirming a change is correct.
golden: build
	@./tests/generate_golden.sh $(BUILD)/pureclip

clean:
	@rm -rf $(BUILD)
