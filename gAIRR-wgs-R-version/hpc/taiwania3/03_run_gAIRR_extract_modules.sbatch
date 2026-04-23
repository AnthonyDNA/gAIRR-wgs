#!/bin/sh
#SBATCH -A CHANGE_ME_SLURM_ACCOUNT    # EDIT_ME: your Taiwania3 project account (e.g. MSTxxxxxx)
#SBATCH -J gAIRR_extract_R
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -o out_extract.log
#SBATCH -e err_extract.log
# Optional: uncomment the next two lines and set your email to get job notifications.
# #SBATCH --mail-user=your.email@example.com
# #SBATCH --mail-type=BEGIN,END,FAIL
# ==============================================================================
# 03_run_gAIRR_extract_modules.sbatch  -  Phase 3 run of gAIRR_extract (R version)
# on a list of BAM / FASTQ samples, using Taiwania3 MODULES (no conda).
#
# Produces per-sample <out>/<sample_id>/on_AIRR_unmapped_R{1,2}.fastq, which you
# can then feed into gAIRR_wgs (Phase 2).
#
# Prerequisites:
#   - 00_setup_modules.sh  (has been run on the login node)
#   - samples.tsv          (see sample_list_template.tsv)
#   - REF_FA               hg38 reference FASTA with BWA index next to it
#                          (files <ref>.fa.bwt, .pac, .sa, .amb, .ann).
#                          If you only have the .fa, run on a login/interactive
#                          node with modules loaded:
#                              bwa index <ref>.fa
#                              samtools faidx <ref>.fa
# ==============================================================================

set -euo pipefail

# ============================== EDIT_ME =======================================
BASE="/staging/biology/$USER"
R_REPO_DIR="$BASE/gAIRR-wgs-R-version"
ORIG_REPO_DIR="$BASE/gAIRR-wgs-main"
VENV_DIR="$BASE/gAIRR-venv"

SAMPLE_LIST="$BASE/gAIRR_R_inputs/samples.tsv"
OUTDIR="$BASE/gAIRR_R_extract_out"
REF_FA="$BASE/reference/hg38_mod.fa"
BED="$ORIG_REPO_DIR/gAIRR_suite/material/extract_materials/gAIRR_allele.bed"
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
for f in "$SAMPLE_LIST" "$BED" "$REF_FA"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: missing file: $f" >&2
        exit 1
    fi
done
if [ ! -f "${REF_FA}.bwt" ]; then
    echo "ERROR: BWA index missing for $REF_FA (expected ${REF_FA}.bwt)." >&2
    echo "       Run:  bwa index $REF_FA && samtools faidx $REF_FA"         >&2
    exit 1
fi

RSCRIPT_MAIN="$R_REPO_DIR/gAIRR_suite_R/gAIRR_extract.R"
if [ ! -f "$RSCRIPT_MAIN" ]; then
    echo "ERROR: R entry script not found: $RSCRIPT_MAIN" >&2
    exit 1
fi
mkdir -p "$OUTDIR"

# --- run ----------------------------------------------------------------------
echo "=== Running gAIRR_extract (R version) ==="
echo "Sample list : $SAMPLE_LIST"
echo "BED         : $BED"
echo "Reference   : $REF_FA"
echo "Output      : $OUTDIR"
echo "Threads     : ${SLURM_CPUS_PER_TASK:-14}"
echo

time Rscript "$RSCRIPT_MAIN" \
    -l "$SAMPLE_LIST" \
    -o "$OUTDIR" \
    -b "$BED" \
    -r "$REF_FA" \
    -t "${SLURM_CPUS_PER_TASK:-14}"

echo
echo "=== Done.  Per-sample outputs under: $OUTDIR/<sample_id>/ ==="
echo "    Each sample folder contains on_AIRR_unmapped_R{1,2}.fastq"
echo "    which you can feed into 02_run_gAIRR_wgs_modules.sbatch."
