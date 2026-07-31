# =============================================================================
# experiments/synthetic_evolution/Makefile
# =============================================================================
#
# Recreates the original `synthetic-c` experiment (../hammock/experiments/
# synthetic-c) against the refactored hammock in this repo. Tests whether:
#   - Mode A best detects Type A (interval insert/delete) changes
#   - Mode B best detects Type B (positional jitter) changes
#   - Mode C best detects Type C (mix of A + B) changes
#
# Bulk artifacts (BEDs, hammock CSVs, logs) live on /vast via symlinks set up
# during scaffolding. Figures and summaries live in the repo.
#
# Usage:
#   make all              # full sweep (≈10 min on one node)
#   make data             # base BED only
#   make evolve           # all evolved snapshot BEDs
#   make hammock          # all hammock CSVs
#   make analyze          # figures + summary CSV
#   make -j 6 hammock     # parallel hammock runs
#   make clean            # rm artifacts in /vast mirror (preserves figures/)
#
# Override RATES / SEED / N_INTERVALS / PRECISION / THREADS on the cmd line:
#   make RATES="0.05" TYPES=A all   # smoke-test with one rate × one type
# =============================================================================

HAMMOCK   ?= /home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock
PYTHON    ?= python3

N_INTERVALS  ?= 100000
SEED         ?= 42
GENERATIONS  ?= 100
STEP         ?= 5
PRECISION    ?= 18
THREADS      ?= 8
# Max start position per chr. Coverage saturation = N*L_avg / (5 * R), so:
#   1e6   → ~47× over-covered (every position covered many times; old default)
#   1e8   → ~0.5% covered (biologically realistic sparse setup; new default)
# Per-genome-range artifacts go to *_g${GENOME_RANGE}/ subdirs so multiple
# runs coexist for side-by-side comparison.
GENOME_RANGE ?= 100000000

RATES        ?= 0.005 0.01 0.05 0.1
TYPES        ?= A B C
# Mode C subB is fixed for the whole experiment (was a sweep; the user
# decided subB only matters in the context of varying expA — so we pin
# subB at one value and sweep expA against it).
SUBB         ?= 0.25
# C_EXPA: Mode C expA sweep — multiplies interval contribution to the sketch.
# expA=0 means "no boost" → CSV name has no _expA suffix.
C_EXPA_POS   ?= 0.50 1.00 2.00 3.00

# ── Paths ───────────────────────────────────────────────────────────────────
CODE     := code
DATA     := data
RESULTS  := results
LOGS     := logs
FIGURES  := figures

# Per-genome-range artifact paths. Tag is "g{GENOME_RANGE}" so artifacts
# from different ranges don't collide.
TAG           := g$(GENOME_RANGE)
BASE_BED      := $(DATA)/synthetic_n$(N_INTERVALS)_$(TAG).bed
BASE_BED_LIST := $(DATA)/synthetic_base_$(TAG).txt
EVOLVED_DIR   := $(DATA)/evolved_$(TAG)
HAMMOCK_DIR   := $(RESULTS)/hammock_$(TAG)
FIGURES_TAG   := $(FIGURES)/$(TAG)
SUMMARY_CSV   := $(RESULTS)/summary_$(TAG).csv
EVOLVE_PREFIX := $(EVOLVED_DIR)/synthetic

# Per-(rate, type): one stamp marking all snapshot BEDs are written, and one
# list-of-paths file fed as bed2 to hammock.
EVOLVE_STAMPS := $(foreach r,$(RATES),$(foreach t,$(TYPES),$(EVOLVED_DIR)/.stamp_r$(r)_$(t)))
LIST_FILES    := $(foreach r,$(RATES),$(foreach t,$(TYPES),$(EVOLVED_DIR)/list_r$(r)_$(t).txt))

# Hammock CSV naming (set by hammock/outprefix.py):
#   {outprefix}_hll_p{P}_jacc{MODE}.csv
#   {outprefix}_hll_p{P}_jaccC[_expA{E}][_B{S}].csv  (suffix omitted if 0/1.0)
# We use outprefix = $(HAMMOCK_DIR)/synth_r{R}_{T} and pin subB at $(SUBB).
PREFIX_FOR    = $(HAMMOCK_DIR)/synth_r$(1)_$(2)
HAMMOCK_AB    := $(foreach r,$(RATES),$(foreach t,$(TYPES),\
                   $(call PREFIX_FOR,$(r),$(t))_hll_p$(PRECISION)_jaccA.csv \
                   $(call PREFIX_FOR,$(r),$(t))_hll_p$(PRECISION)_jaccB.csv))
# Mode C, expA=0 (no boost), subB=$(SUBB) → name has _B{SUBB} but no _expA.
HAMMOCK_C_DEF := $(foreach r,$(RATES),$(foreach t,$(TYPES),\
                   $(call PREFIX_FOR,$(r),$(t))_hll_p$(PRECISION)_jaccC_B$(SUBB).csv))
# Mode C, expA > 0, subB=$(SUBB) → name has both _expA{E} and _B{SUBB}.
HAMMOCK_C_EXP := $(foreach r,$(RATES),$(foreach t,$(TYPES),$(foreach e,$(C_EXPA_POS),\
                   $(call PREFIX_FOR,$(r),$(t))_hll_p$(PRECISION)_jaccC_expA$(e)_B$(SUBB).csv)))
HAMMOCK_CSVS  := $(HAMMOCK_AB) $(HAMMOCK_C_DEF) $(HAMMOCK_C_EXP)

FIGURE_TARGETS  := \
  $(FIGURES_TAG)/jaccard_by_generation.png \
  $(FIGURES_TAG)/sensitivity_heatmap.png \
  $(FIGURES_TAG)/mode_c_expA_sweep.png \
  $(FIGURES_TAG)/mode_c_expA_interpolation.png \
  $(foreach r,$(RATES),$(FIGURES_TAG)/mode_per_type_r$(r).png)
SUMMARY_TARGETS := $(SUMMARY_CSV)

# Allow stems containing dots and one underscore (rate_type → e.g. "0.05_A").
.SECONDEXPANSION:

# ── Top-level targets ───────────────────────────────────────────────────────
.PHONY: all data evolve hammock analyze clean help list-targets

all: analyze

data: $(BASE_BED) $(BASE_BED_LIST)

evolve: $(EVOLVE_STAMPS) $(LIST_FILES)

hammock: $(HAMMOCK_CSVS)

analyze: $(FIGURE_TARGETS) $(SUMMARY_TARGETS)

help:
	@echo "Targets: all, data, evolve, hammock, analyze, clean"
	@echo "Vars:    RATES=\"$(RATES)\" TYPES=\"$(TYPES)\" SUBB=$(SUBB) C_EXPA_POS=\"$(C_EXPA_POS)\""
	@echo "         N_INTERVALS=$(N_INTERVALS) PRECISION=$(PRECISION) THREADS=$(THREADS)"
	@echo "         GENOME_RANGE=$(GENOME_RANGE) → tag '$(TAG)'"
	@echo "Artifacts:"
	@echo "         data/synthetic_n${N_INTERVALS}_$(TAG).bed (base)"
	@echo "         $(EVOLVED_DIR)/  (snapshots)"
	@echo "         $(HAMMOCK_DIR)/  (hammock CSVs)"
	@echo "         $(FIGURES_TAG)/  (figures)"
	@echo "         $(SUMMARY_CSV)   (summary)"

# ── Step 1: base BED ────────────────────────────────────────────────────────
$(BASE_BED): $(CODE)/createBED.py
	$(PYTHON) $(CODE)/createBED.py \
	    --num_lines $(N_INTERVALS) \
	    --output_file $@ \
	    --genome_range $(GENOME_RANGE) \
	    --seed $(SEED)

$(BASE_BED_LIST): $(BASE_BED)
	realpath $(BASE_BED) > $@

# ── Step 2: evolved snapshots (one evolve.py call per rate × type) ──────────
# Stem: {rate}_{type}, e.g. "0.05_A". evolve.py writes 20 BEDs at G=100 step=5.
$(EVOLVED_DIR)/.stamp_r%: $(BASE_BED) $(CODE)/evolve.py
	@mkdir -p $(EVOLVED_DIR) $(LOGS)
	@rate=$$(echo "$*" | sed -E 's/_([ABC])$$//'); \
	type=$$(echo "$*" | sed -E 's/.*_//'); \
	echo "[evolve $(TAG)] rate=$$rate type=$$type → $(EVOLVED_DIR)/synthetic_g$(GENERATIONS)_r$$rate\_$$type\_*.bed"; \
	$(PYTHON) $(CODE)/evolve.py $(BASE_BED) $(EVOLVE_PREFIX) \
	    --generations $(GENERATIONS) \
	    --prob $$rate \
	    --type $$type \
	    --step $(STEP) \
	    --seed $(SEED) \
	    --genome_range $(GENOME_RANGE) \
	    > $(LOGS)/evolve_$(TAG)_r$$rate\_$$type.log 2>&1
	touch $@

# List file (bed2 for hammock), one per (rate, type).
$(EVOLVED_DIR)/list_r%.txt: $(EVOLVED_DIR)/.stamp_r$$*
	@rate=$$(echo "$*" | sed -E 's/_([ABC])$$//'); \
	type=$$(echo "$*" | sed -E 's/.*_//'); \
	ls -1 $(EVOLVED_DIR)/synthetic_g$(GENERATIONS)_r$$rate\_$$type\_*.bed \
	  | xargs realpath > $@

# ── Step 3: hammock runs ────────────────────────────────────────────────────
# Mode A
$(HAMMOCK_DIR)/synth_r%_hll_p$(PRECISION)_jaccA.csv: \
        $(BASE_BED_LIST) \
        $(EVOLVED_DIR)/list_r$$*.txt
	@mkdir -p $(HAMMOCK_DIR) $(LOGS)
	@rate=$$(echo "$*" | sed -E 's/_([ABC])$$//'); \
	type=$$(echo "$*" | sed -E 's/.*_//'); \
	echo "[hammock A] rate=$$rate type=$$type"; \
	$(HAMMOCK) $(BASE_BED_LIST) $(EVOLVED_DIR)/list_r$*.txt \
	    --mode A --precision $(PRECISION) \
	    --threads $(THREADS) \
	    --outprefix $(HAMMOCK_DIR)/synth_r$$rate\_$$type \
	    > $(LOGS)/hammock_r$$rate\_$$type\_modeA.log 2>&1

# Mode B
$(HAMMOCK_DIR)/synth_r%_hll_p$(PRECISION)_jaccB.csv: \
        $(BASE_BED_LIST) \
        $(EVOLVED_DIR)/list_r$$*.txt
	@mkdir -p $(HAMMOCK_DIR) $(LOGS)
	@rate=$$(echo "$*" | sed -E 's/_([ABC])$$//'); \
	type=$$(echo "$*" | sed -E 's/.*_//'); \
	echo "[hammock B] rate=$$rate type=$$type"; \
	$(HAMMOCK) $(BASE_BED_LIST) $(EVOLVED_DIR)/list_r$*.txt \
	    --mode B --precision $(PRECISION) \
	    --threads $(THREADS) \
	    --outprefix $(HAMMOCK_DIR)/synth_r$$rate\_$$type \
	    > $(LOGS)/hammock_r$$rate\_$$type\_modeB.log 2>&1

# Mode C, expA=0 (no boost), subB=$(SUBB). CSV name has _B{SUBB} but no _expA.
# Static pattern: one CSV per (rate, type).
$(HAMMOCK_C_DEF): $(BASE_BED_LIST)
	@mkdir -p $(HAMMOCK_DIR) $(LOGS)
	@bn=$$(basename $@ .csv); \
	rate=$$(echo "$$bn" | sed -E 's/^synth_r([0-9.]+)_[ABC]_.*/\1/'); \
	type=$$(echo "$$bn" | sed -E 's/^synth_r[0-9.]+_([ABC])_.*/\1/'); \
	list=$(EVOLVED_DIR)/list_r$$rate\_$$type.txt; \
	$(MAKE) --no-print-directory $$list; \
	echo "[hammock C subB=$(SUBB) expA=0] rate=$$rate type=$$type"; \
	$(HAMMOCK) $(BASE_BED_LIST) $$list \
	    --mode C --precision $(PRECISION) --subB $(SUBB) \
	    --threads $(THREADS) \
	    --outprefix $(HAMMOCK_DIR)/synth_r$$rate\_$$type \
	    > $(LOGS)/hammock_r$$rate\_$$type\_modeC_B$(SUBB).log 2>&1

# Mode C, expA > 0, subB=$(SUBB). Static pattern, one CSV per (rate, type, expA).
# Filename: synth_r{R}_{T}_hll_p{P}_jaccC_expA{E}_B{SUBB}.csv
$(HAMMOCK_C_EXP): $(BASE_BED_LIST)
	@mkdir -p $(HAMMOCK_DIR) $(LOGS)
	@bn=$$(basename $@ .csv); \
	rate=$$(echo "$$bn" | sed -E 's/^synth_r([0-9.]+)_[ABC]_.*/\1/'); \
	type=$$(echo "$$bn" | sed -E 's/^synth_r[0-9.]+_([ABC])_.*/\1/'); \
	expa=$$(echo "$$bn" | sed -E 's/.*_expA([0-9.]+)_B.*/\1/'); \
	list=$(EVOLVED_DIR)/list_r$$rate\_$$type.txt; \
	$(MAKE) --no-print-directory $$list; \
	echo "[hammock C expA=$$expa subB=$(SUBB)] rate=$$rate type=$$type"; \
	$(HAMMOCK) $(BASE_BED_LIST) $$list \
	    --mode C --precision $(PRECISION) --expA $$expa --subB $(SUBB) \
	    --threads $(THREADS) \
	    --outprefix $(HAMMOCK_DIR)/synth_r$$rate\_$$type \
	    > $(LOGS)/hammock_r$$rate\_$$type\_modeC_expA$$expa\_B$(SUBB).log 2>&1

# ── Step 4: analysis (figures + summary CSV) ────────────────────────────────
$(FIGURE_TARGETS) $(SUMMARY_TARGETS) &: $(HAMMOCK_CSVS) $(CODE)/analyze.R
	@mkdir -p $(FIGURES_TAG) $(LOGS) $(dir $(SUMMARY_TARGETS))
	@command -v Rscript >/dev/null 2>&1 || { \
	  echo "Rscript not on PATH. On Rockfish: ml r/4.3.0"; exit 1; }
	Rscript $(CODE)/analyze.R \
	    --hammock-dir $(HAMMOCK_DIR) \
	    --out-figures $(FIGURES_TAG) \
	    --out-summary $(SUMMARY_TARGETS) \
	    > $(LOGS)/analyze_$(TAG).log 2>&1

# ── Cleanup ─────────────────────────────────────────────────────────────────
# Wipes the per-tag evolved snapshots, hammock CSVs, and logs for the
# CURRENT GENOME_RANGE. Other ranges (data/evolved_g*/, results/hammock_g*/)
# are left alone. Base BED and figures/ also preserved.
clean:
	rm -rf $(EVOLVED_DIR)/* $(HAMMOCK_DIR)/* $(LOGS)/*_$(TAG)_*
	@echo "Removed $(TAG) evolved snapshots, hammock CSVs, and tagged logs."
	@echo "Kept base BED ($(BASE_BED)), figures/, and other genome ranges."

list-targets:
	@echo "EVOLVE_STAMPS:"; for f in $(EVOLVE_STAMPS); do echo "  $$f"; done
	@echo "LIST_FILES:";    for f in $(LIST_FILES);    do echo "  $$f"; done
	@echo "HAMMOCK_CSVS:";  for f in $(HAMMOCK_CSVS);  do echo "  $$f"; done
