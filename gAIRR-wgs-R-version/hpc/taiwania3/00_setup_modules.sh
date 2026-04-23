#!/bin/bash
# ==============================================================================
# 00_setup_modules.sh  -  one-time setup for the R version of gAIRR-wgs on
# Taiwania3 using the MODULE system (no conda).
#
# Run this ONCE on the LOGIN NODE (not as an sbatch job).
#
# Assumes you've already placed the R-version repo at:
#     /staging/biology/$USER/gAIRR-wgs-R-version
# and test FASTQs at:
#     /staging/biology/$USER/gAIRR-wgs-R-version/test_sample/HG002_part_R{1,2}.fastq.gz
#
# Steps:
#   1. module-load the toolchain (R, BWA, Samtools, Python, SPAdes).
#   2. Clone the ORIGINAL gAIRR-wgs repo next to the R version (if missing).
#      We need it because the R port reuses its scripts/ and material/.
#   3. Create a Python virtualenv for the helper-script deps (pysam, pyfastx,
#      numpy, pandas, matplotlib).
#   4. Symlink scripts/ and material/ into the R tree via clone_script.sh.
#   5. Sanity check: show which tools resolve to which paths and run
#      `gAIRR_wgs --help` via Rscript to confirm the port parses.
# ==============================================================================

set -euo pipefail

# ----- adjust these if you chose different paths ------------------------------
BASE="/staging/biology/$USER"
R_REPO_DIR="$BASE/gAIRR-wgs-R-version"
ORIG_REPO_DIR="$BASE/gAIRR-wgs-main"
VENV_DIR="$BASE/gAIRR-venv"
# ------------------------------------------------------------------------------

echo "==> Loading Taiwania3 modules..."
module load biology
module load R/4.4.1
module load BWA/0.7.17
module load Samtools/1.15.1
module load python/3.12.2
module load SPAdes/3.15.5

echo "==> Modules loaded. Tool paths:"
for t in Rscript bwa samtools python3 spades.py; do
    printf "  %-12s %s\n" "$t" "$(command -v "$t" 2>/dev/null || echo MISSING)"
done

# ----- 1. clone original repo -------------------------------------------------
if [ ! -d "$ORIG_REPO_DIR/gAIRR_suite" ]; then
    echo
    echo "==> Cloning original gAIRR-wgs into $ORIG_REPO_DIR ..."
    git clone https://github.com/AnthonyDNA/gAIRR-wgs.git "$ORIG_REPO_DIR"
else
    echo
    echo "==> Original gAIRR-wgs already present at $ORIG_REPO_DIR (skipping clone)"
fi

# ----- 2. Python venv for helper-script deps ----------------------------------
if [ ! -d "$VENV_DIR" ]; then
    echo
    echo "==> Creating Python venv at $VENV_DIR ..."
    python3 -m venv "$VENV_DIR"
fi

# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

echo "==> Installing Python deps (this is a one-off; takes a minute)..."
pip install --upgrade pip > /dev/null
pip install \
    "pysam>=0.22,<0.23" \
    "pyfastx>=2.0,<3" \
    "numpy" \
    "pandas" \
    "matplotlib"

echo "==> pip packages installed:"
pip list 2>/dev/null | grep -Ei '^(pysam|pyfastx|numpy|pandas|matplotlib)\b'

# ----- 3. symlink scripts/ and material/ via clone_script.sh ------------------
echo
echo "==> Wiring up scripts/ and material/ symlinks..."
cd "$R_REPO_DIR"
GAIRR_ORIG="$ORIG_REPO_DIR" bash clone_script.sh

# ----- 4. sanity check --------------------------------------------------------
echo
echo "======================== sanity check ========================"
echo "Rscript --help on the port:"
Rscript "$R_REPO_DIR/gAIRR_suite_R/gAIRR_wgs.R" --help | head -6
echo
echo "Helper-script import test (confirms venv can see pysam/pyfastx):"
python3 -c "import pysam, pyfastx, numpy, pandas, matplotlib; print('  OK: all helper imports work')"
echo "=============================================================="
echo
echo "Setup OK."
echo "Next step:  sbatch 01_smoke_test_modules.sbatch"
