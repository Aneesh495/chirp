# chirp - High-Resolution FMCW 4D Imaging Radar & SAR Engine
# Build automation for MEX acceleration kernels, testing, and benchmarks

OCTAVE ?= octave
OCTAVE_CLI ?= octave-cli --no-gui
MKOCTFILE ?= mkoctfile
MEX_DIR = src/mex
MEX_OUT = bin/mex

.PHONY: all mex test audit bench clean demo help

all: mex test

help:
	@echo "chirp - High-Resolution FMCW 4D Imaging Radar Perception Engine"
	@echo ""
	@echo "Targets:"
	@echo "  mex      Compile native C MEX acceleration kernels"
	@echo "  test     Execute full unit and integration test suite"
	@echo "  audit    Execute independent mathematical invariants and negative controls"
	@echo "  bench    Run computational performance benchmarks"
	@echo "  demo     Run 4D imaging radar perception demonstration"
	@echo "  clean    Remove compiled binaries and intermediate build artifacts"

$(MEX_OUT):
	mkdir -p $(MEX_OUT)

mex: $(MEX_OUT)
	@echo "Compiling C MEX acceleration kernels..."
	$(MKOCTFILE) --mex $(MEX_DIR)/chirp_cfar_2d_mex.c -o $(MEX_OUT)/chirp_cfar_2d_mex.mex
	$(MKOCTFILE) --mex $(MEX_DIR)/chirp_range_fft_mex.c -o $(MEX_OUT)/chirp_range_fft_mex.mex
	$(MKOCTFILE) --mex $(MEX_DIR)/chirp_backprojection_mex.c -o $(MEX_OUT)/chirp_backprojection_mex.mex
	@echo "MEX acceleration kernels compiled successfully."

test:
	@echo "Running chirp automated test suite..."
	$(OCTAVE_CLI) --eval "addpath('src'); addpath('tests'); addpath('bin/mex'); run_all_tests;"

audit:
	@echo "Running chirp mathematical invariant audit and negative controls..."
	$(OCTAVE_CLI) --eval "addpath('src'); addpath('tests'); addpath('bin/mex'); test_negative_controls; test_mathematical_invariants; test_no_em_dashes;"

bench:
	@echo "Running chirp computational benchmarks..."
	$(OCTAVE_CLI) --eval "addpath('src'); addpath('benchmarks'); addpath('bin/mex'); run_benchmarks;"

demo:
	@echo "Running chirp 4D automotive perception demo..."
	$(OCTAVE_CLI) --eval "addpath('src'); addpath('examples'); addpath('bin/mex'); demo_4d_automotive_radar;"

clean:
	rm -rf $(MEX_OUT) output/ test_output/ build/
	find . -type f -name "*.mex" -delete
	find . -type f -name "*.o" -delete
	find . -type f -name "*.asv" -delete
	find . -type f -name "octave-core" -delete
