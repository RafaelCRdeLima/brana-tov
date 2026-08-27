CXX ?= g++
CXXFLAGS ?= -std=c++17 -Wall -Wextra -pedantic -O2 -fopenmp
CPPFLAGS ?= -Iincludes
LDFLAGS ?= -fopenmp

TARGET := build/brana_tov
SOURCES := src/main.cpp src/eos.cpp src/mcmc.cpp src/pchip.cpp src/tov.cpp
STEPS ?= 800
WALKERS ?= 32
ENSEMBLES ?= 8
THREADS ?= 8
PROFILE ?= default
EOS ?= data/SLy.txt
OUT ?= output/data
LAMBDA ?= 1e3
ALPHA ?= -0.25
W ?= 0.3333333333
MAMBA_ENV ?= bwmcmc
PYTHON := micromamba run -n $(MAMBA_ENV) python

.PHONY: all clean sequence mcmc plots modified-tov corner campaign-sequences campaign-smoke campaign-full campaign-full-remaining campaign-summary campaign-paper

all: $(TARGET)

GIT_COMMIT := $(shell git rev-parse HEAD 2>/dev/null || echo unknown)
GIT_DIRTY  := $(shell test -z "`git status --porcelain 2>/dev/null`" && echo clean || echo dirty)
PROVDEFS   := -DBRANA_GIT_COMMIT='"$(GIT_COMMIT)"' -DBRANA_GIT_DIRTY='"$(GIT_DIRTY)"' \
              -DBRANA_BUILD_FLAGS='"$(CXXFLAGS)"'

$(TARGET): $(SOURCES)
	mkdir -p build
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(PROVDEFS) $(SOURCES) $(LDFLAGS) -o $(TARGET)

sequence: $(TARGET)
	./$(TARGET) sequence $(EOS) $(OUT)

mcmc: $(TARGET)
	./$(TARGET) mcmc $(STEPS) $(WALKERS) $(ENSEMBLES) $(THREADS) $(PROFILE) $(EOS) $(OUT)

campaign-sequences: $(TARGET)
	./$(TARGET) sequence data/SLy.txt output/campaign/sequences/SLy
	./$(TARGET) sequence data/eos/tables/BSK24.txt output/campaign/sequences/BSK24
	./$(TARGET) sequence data/eos/tables/DD2.txt output/campaign/sequences/DD2
	./$(TARGET) sequence data/eos/tables/NL3.txt output/campaign/sequences/NL3
	./$(TARGET) sequence data/eos/tables/QHC19A.txt output/campaign/sequences/QHC19A
	./$(TARGET) sequence data/eos/tables/QHC21D.txt output/campaign/sequences/QHC21D
	./$(TARGET) sequence data/eos/tables/SPGM4.txt output/campaign/sequences/SPGM4

campaign-smoke: $(TARGET)
	$(MAKE) mcmc STEPS=80 WALKERS=16 ENSEMBLES=4 THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/SLy.txt OUT=output/campaign/smoke/SLy
	$(MAKE) mcmc STEPS=80 WALKERS=16 ENSEMBLES=4 THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/BSK24.txt OUT=output/campaign/smoke/BSK24
	$(MAKE) mcmc STEPS=80 WALKERS=16 ENSEMBLES=4 THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/DD2.txt OUT=output/campaign/smoke/DD2
	$(MAKE) mcmc STEPS=80 WALKERS=16 ENSEMBLES=4 THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/NL3.txt OUT=output/campaign/smoke/NL3
	$(MAKE) mcmc STEPS=80 WALKERS=16 ENSEMBLES=4 THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/QHC19A.txt OUT=output/campaign/smoke/QHC19A
	$(MAKE) mcmc STEPS=80 WALKERS=16 ENSEMBLES=4 THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/QHC21D.txt OUT=output/campaign/smoke/QHC21D
	$(MAKE) mcmc STEPS=80 WALKERS=16 ENSEMBLES=4 THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/SPGM4.txt OUT=output/campaign/smoke/SPGM4

campaign-full: $(TARGET)
	$(MAKE) mcmc STEPS=$(STEPS) WALKERS=$(WALKERS) ENSEMBLES=$(ENSEMBLES) THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/SLy.txt OUT=output/campaign/full/SLy
	$(MAKE) campaign-full-remaining STEPS=$(STEPS) WALKERS=$(WALKERS) ENSEMBLES=$(ENSEMBLES) THREADS=$(THREADS)

campaign-full-remaining: $(TARGET)
	$(MAKE) mcmc STEPS=$(STEPS) WALKERS=$(WALKERS) ENSEMBLES=$(ENSEMBLES) THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/BSK24.txt OUT=output/campaign/full/BSK24
	$(MAKE) mcmc STEPS=$(STEPS) WALKERS=$(WALKERS) ENSEMBLES=$(ENSEMBLES) THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/DD2.txt OUT=output/campaign/full/DD2
	$(MAKE) mcmc STEPS=$(STEPS) WALKERS=$(WALKERS) ENSEMBLES=$(ENSEMBLES) THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/NL3.txt OUT=output/campaign/full/NL3
	$(MAKE) mcmc STEPS=$(STEPS) WALKERS=$(WALKERS) ENSEMBLES=$(ENSEMBLES) THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/QHC19A.txt OUT=output/campaign/full/QHC19A
	$(MAKE) mcmc STEPS=$(STEPS) WALKERS=$(WALKERS) ENSEMBLES=$(ENSEMBLES) THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/QHC21D.txt OUT=output/campaign/full/QHC21D
	$(MAKE) mcmc STEPS=$(STEPS) WALKERS=$(WALKERS) ENSEMBLES=$(ENSEMBLES) THREADS=$(THREADS) PROFILE=perturbative-gr EOS=data/eos/tables/SPGM4.txt OUT=output/campaign/full/SPGM4

campaign-summary:
	python3 scripts/summarize_eos_campaign.py

campaign-paper:
	python3 scripts/check_eos_campaign.py
	python3 scripts/summarize_eos_campaign.py --mcmc output/campaign/full --output output/campaign/full/summary.csv
	$(PYTHON) scripts/build_eos_paper_assets.py --summary output/campaign/full/summary.csv --sequences output/campaign/sequences --output-dir paper2 --campaign-label production --diagnostics output/campaign/full/diagnostics.csv
	$(PYTHON) scripts/plot_multi_eos_mass_radius.py --summary output/campaign/full/summary.csv --sequences output/campaign/sequences --median-dir output/campaign/full/median_sequences --output paper2/multi_eos_mass_radius.png

plots:
	$(PYTHON) scripts/plot_sequences.py
	$(PYTHON) scripts/plot_mcmc.py

modified-tov: $(TARGET)
	$(PYTHON) scripts/plot_modified_tov.py --lambda $(LAMBDA) --alpha $(ALPHA) --w $(W) --no-build

corner:
	$(PYTHON) scripts/plot_corner.py

clean:
	rm -rf build
