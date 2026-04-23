#!/bin/sh
#SBATCH -A CHANGE_ME_SLURM_ACCOUNT    # EDIT_ME: your Taiwania3 project account (e.g. MSTxxxxxx)
#SBATCH -J gAIRRwgs_R
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -o out.log
#SBATCH -e err.log
# Optional: uncomment the next two lines and set your email to get job notifications.
# #SBATCH --mail-user=your.email@example.com
# #SBATCH --mail-type=BEGIN,END,FAIL
# ==============================================================================
# 02_run_gAIRR_wgs_modules.sbatch  -  Phase 2 run of gAIRR_wgs (R version) on a
# real WGS sample, using Taiwania3 MODULES (no conda).
#
# Prerequisite: bash 00_setup_modules.sh  (run once on the login node)
#
# For a 30x WGS sample, expect ~1-3 hours for all six loci without flanking,
# considerably longer with --flanking (SPAdes assembly).
# ==============================================================================

set -euo pipefail

# ============================== EDIT_ME =======================================
BASE="/staging/biology/$USER"
R_REPO_DIR="$BASE/gAIRR-wgs-R-version"
VENV_DIR="$BASE/gAIRR-venv"

SAMPLE_ID="HG002"
READ1="/staging/biology/$USER/wgs_fastq/HG002.R1.fastq.gz"
READ2="/staging/biology/$USER/wgs_fastq/HG002.R2.fastq.gz"
WORK_DIR="/staging/biology/$USER/gAIRR_R_out/target_call"
LOCI="TRV TRJ TRD IGV IGJ IGD"      # subset: e.g. LOCI="TRV TRJ TRD"
FLANKING=""                          # set FLANKING="--flanking" to enable
KEEP=""                              # set KEEP="--keep" to preserve intermediate .sam
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

# --- run ----------------------------------------------------------------------
echo "=== Running gAIRR_wgs (R version) ==="
echo "Sample     : $SAMPLE_ID"
echo "R1         : $READ1"
echo "R2         : $READ2"
echo "Loci       : $LOCI"
echo "Flanking   : ${FLANKING:-no}"
echo "Keep SAMs  : ${KEEP:-no}"
echo "Threads    : ${SLURM_CPUS_PER_TASK:-14}"
echo "Output dir : $WORK_DIR"
echo

time Rscript "$RSCRIPT_MAIN" \
    -wd  "$WORK_DIR" \
    -id  "$SAMPLE_ID" \
    -rd1 "$READ1" \
    -rd2 "$READ2" \
    -lc  $LOCI \
    -t   "${SLURM_CPUS_PER_TASK:-14}" \
    $FLANKING \
    $KEEP

echo
echo "=== Done.  Results under: $WORK_DIR/$SAMPLE_ID/ ==="
