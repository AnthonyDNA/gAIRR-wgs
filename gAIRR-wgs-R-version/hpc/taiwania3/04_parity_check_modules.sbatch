#!/bin/bash
#SBATCH -A CHANGE_ME_SLURM_ACCOUNT    # EDIT_ME: your Taiwania3 project account (e.g. MSTxxxxxx)
#SBATCH -J gAIRRwgs_parity
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -o out_parity.log
#SBATCH -e err_parity.log
# Optional: uncomment the next two lines and set your email to get job notifications.
# #SBATCH --mail-user=your.email@example.com
# #SBATCH --mail-type=BEGIN,END,FAIL
# ==============================================================================
# 04_parity_check_modules.sbatch  -  Phase 4 parity check using Taiwania3
# MODULES (no conda).
#
# Runs the ORIGINAL Python version and the R version on the same FASTQ pair,
# into parallel output trees, then diffs them.  The R port delegates all
# sequence work to the same underlying scripts, so diffs here should be zero.
#
# Prerequisite: bash 00_setup_modules.sh  (run once on the login node)
# ==============================================================================

set -euo pipefail

# ============================== EDIT_ME =======================================
BASE="/staging/biology/$USER"
R_REPO_DIR="$BASE/gAIRR-wgs-R-version"
ORIG_REPO_DIR="$BASE/gAIRR-wgs-main"
VENV_DIR="$BASE/gAIRR-venv"

SAMPLE_ID="HG002_parity"
READ1="$R_REPO_DIR/test_sample/HG002_part_R1.fastq.gz"
READ2="$R_REPO_DIR/test_sample/HG002_part_R2.fastq.gz"
OUT_BASE="$BASE/gAIRR_parity"
LOCI="TRV TRJ TRD IGV IGJ IGD"                      
# =============================== end EDIT_ME ==================================

# --- load modules -------------------------------------------------------------
module load biology
module load R/4.4.1
module load BWA/0.7.17
module load Samtools/1.15.1
module load python/3.12.2
module load SPAdes/3.15.5

# --- activate python venv -----------------------------------------------------
if [ ! -d "$VENV_DIR" ]; then
    echo "ERROR: Python venv missing at $VENV_DIR" >&2
    echo "       Run 00_setup_modules.sh on the login node first." >&2
    exit 1
fi
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

# --- validate -----------------------------------------------------------------
if [ ! -d "$ORIG_REPO_DIR/gAIRR_suite" ]; then
    echo "ERROR: original gAIRR-wgs not found at $ORIG_REPO_DIR" >&2
    exit 1
fi
for f in "$READ1" "$READ2"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: FASTQ not found: $f" >&2
        exit 1
    fi
done

PY_OUT="$OUT_BASE/py"
R_OUT="$OUT_BASE/r"
mkdir -p "$PY_OUT" "$R_OUT"

T="${SLURM_CPUS_PER_TASK:-14}"

# --- 1. Python version --------------------------------------------------------
echo "=============================================================="
echo " Running PYTHON version (original)"
echo "=============================================================="
time python3 "$ORIG_REPO_DIR/gAIRR_suite/gAIRR_wgs.py" \
    -wd  "$PY_OUT" \
    -id  "$SAMPLE_ID" \
    -rd1 "$READ1" \
    -rd2 "$READ2" \
    -lc  $LOCI \
    -t   "$T"

# --- 2. R version -------------------------------------------------------------
echo "=============================================================="
echo " Running R version"
echo "=============================================================="
time Rscript "$R_REPO_DIR/gAIRR_suite_R/gAIRR_wgs.R" \
    -wd  "$R_OUT" \
    -id  "$SAMPLE_ID" \
    -rd1 "$READ1" \
    -rd2 "$READ2" \
    -lc  $LOCI \
    -t   "$T"

# --- 3. diff ------------------------------------------------------------------
echo "=============================================================="
echo " Diffing outputs"
echo "=============================================================="

PY_DIR="$PY_OUT/$SAMPLE_ID"
R_DIR="$R_OUT/$SAMPLE_ID"

fail=0

py_list=$(cd "$PY_DIR" && find . -type f | LC_ALL=C sort)
r_list=$( cd "$R_DIR"  && find . -type f | LC_ALL=C sort)
if [ "$py_list" != "$r_list" ]; then
    echo "DIFF: file lists differ ----"
    diff <(echo "$py_list") <(echo "$r_list") | head -80
    fail=$((fail+1))
else
    echo "OK: same set of output files ($(echo "$py_list" | wc -l) files)"
fi

while IFS= read -r rel; do
    a="$PY_DIR/$rel"
    b="$R_DIR/$rel"
    if [ ! -f "$a" ] || [ ! -f "$b" ]; then continue; fi
    ma=$(md5sum "$a" | awk '{print $1}')
    mb=$(md5sum "$b" | awk '{print $1}')
    if [ "$ma" = "$mb" ]; then
        printf "  match  %s\n" "$rel"
    else
        printf "  DIFF   %s  (py=%s  r=%s)\n" "$rel" "$ma" "$mb"
        fail=$((fail+1))
    fi
done < <(cd "$PY_DIR" && find . -type f \( -name '*.fasta' -o -name '*.rpt' -o -name 'bwa_alleles_all.txt' \) | LC_ALL=C sort)

if [ "$fail" -eq 0 ]; then
    echo
    echo "=============================================================="
    echo " PARITY CHECK PASSED - R and Python outputs match."
    echo "=============================================================="
else
    echo
    echo "=============================================================="
    echo " PARITY CHECK FOUND $fail DIFFERENCE(S) - investigate above." >&2
    echo "=============================================================="
    exit 2
fi
