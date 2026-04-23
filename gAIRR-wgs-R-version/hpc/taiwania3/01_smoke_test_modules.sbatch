#!/bin/sh
#SBATCH -A CHANGE_ME_SLURM_ACCOUNT    # EDIT_ME: your Taiwania3 project account (e.g. MSTxxxxxx)
#SBATCH -J gAIRRwgs_R_smoke
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -o out_smoke.log
#SBATCH -e err_smoke.log
# Optional: uncomment the next two lines and set your email to get job notifications.
# #SBATCH --mail-user=your.email@example.com
# #SBATCH --mail-type=BEGIN,END,FAIL
# ==============================================================================
# 01_smoke_test_modules.sbatch  -  Phase 1 smoke test of the R version using the
# Taiwania3 MODULE system (no conda).
#
# Runs the R port on the test FASTQs you staged at:
#   /staging/biology/$USER/gAIRR-wgs-R-version/test_sample/HG002_part_R1.fastq.gz
#   /staging/biology/$USER/gAIRR-wgs-R-version/test_sample/HG002_part_R2.fastq.gz
#
# Uses only TRV (fastest locus) so it finishes in ~5 minutes.  A pass here means
# the R wrapper parses its args, resolves all inner scripts/tools, and completes
# at least one full locus end-to-end.
#
# Prerequisite:
#   bash 00_setup_modules.sh   (run once on the login node)
# ==============================================================================

set -euo pipefail

# ============================== EDIT_ME =======================================
BASE="/staging/biology/$USER"
R_REPO_DIR="$BASE/gAIRR-wgs-R-version"
VENV_DIR="$BASE/gAIRR-venv"

READ1="$R_REPO_DIR/test_sample/HG002_part_R1.fastq.gz"
READ2="$R_REPO_DIR/test_sample/HG002_part_R2.fastq.gz"

WORK_DIR="$R_REPO_DIR/hpc/taiwania3/smoke_out"
SAMPLE_ID="HG002_smoke"
LOCI="TRV TRJ TRD IGV IGJ IGD"                    # all six loci; set to e.g. "TRV" for a ~1-minute smoke test
# =============================== end EDIT_ME ==================================

# --- load modules -------------------------------------------------------------
echo "==> Loading Taiwania3 modules..."
module load biology
module load R/4.4.1
module load BWA/0.7.17
module load Samtools/1.15.1
module load python/3.12.2
module load SPAdes/3.15.5

# --- activate python venv for helper deps (pysam / pyfastx / ...) -------------
if [ ! -d "$VENV_DIR" ]; then
    echo "ERROR: Python venv missing at $VENV_DIR" >&2
    echo "       Run 00_setup_modules.sh on the login node first." >&2
    exit 1
fi
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

# --- validate -----------------------------------------------------------------
for f in "$READ1" "$READ2"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: FASTQ not found: $f" >&2
        exit 1
    fi
done

RSCRIPT_MAIN="$R_REPO_DIR/gAIRR_suite_R/gAIRR_wgs.R"
if [ ! -f "$RSCRIPT_MAIN" ]; then
    echo "ERROR: R entry script not found: $RSCRIPT_MAIN" >&2
    exit 1
fi

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# --- sanity print -------------------------------------------------------------
echo
echo "=== R version sanity check ==="
for t in Rscript bwa samtools python3 spades.py; do
    printf "  %-12s %s\n" "$t" "$(command -v "$t" 2>/dev/null || echo MISSING)"
done
Rscript --version
echo

# --- run ----------------------------------------------------------------------
echo "=== Running gAIRR_wgs (R version) on test_sample ==="
echo "Sample     : $SAMPLE_ID"
echo "R1         : $READ1"
echo "R2         : $READ2"
echo "Loci       : $LOCI"
echo "Output dir : $WORK_DIR"
echo "Threads    : ${SLURM_CPUS_PER_TASK:-14}"
echo

time Rscript "$RSCRIPT_MAIN" \
    -wd  "$WORK_DIR" \
    -id  "$SAMPLE_ID" \
    -rd1 "$READ1" \
    -rd2 "$READ2" \
    -lc  $LOCI \
    -t   "${SLURM_CPUS_PER_TASK:-14}"

# --- verify outputs -----------------------------------------------------------
echo
echo "=== Verifying outputs ==="
OUT="$WORK_DIR/$SAMPLE_ID"
if [ ! -d "$OUT" ]; then
    echo "FAIL: expected output directory $OUT not found" >&2
    exit 1
fi

EXPECT_FILES=(
    "$OUT/${SAMPLE_ID}_TCRV_novel/TCRV_with_novel.fasta"
    "$OUT/${SAMPLE_ID}_TCRV"
)
missing=0
for f in "${EXPECT_FILES[@]}"; do
    if [ ! -e "$f" ]; then
        echo "  MISSING: $f" >&2
        missing=$((missing+1))
    else
        echo "  OK: $f"
    fi
done
if [ "$missing" -gt 0 ]; then
    echo "FAIL: $missing expected output(s) missing" >&2
    exit 1
fi

echo
echo "=== Smoke test PASSED ==="
echo "Full output tree under: $OUT"
