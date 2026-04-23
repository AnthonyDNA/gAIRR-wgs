# gAIRR-wgs (R version)

An R port of the user-facing entry points of [gAIRR-wgs](https://github.com/AnthonyDNA/gAIRR-wgs). The original project is a bioinformatics pipeline that genotypes germline **T-cell receptor (TCR)** and **immunoglobulin (IG/BCR)** genes from standard whole-genome sequencing data (paired-end 150 bp, ≥30× coverage). This R version exposes the same command-line interface in R and produces identical outputs across all six gAIRR loci (TCRV, TCRJ, TCRD_plusHep, BCRV, BCRJ, BCRD_plusHep).

## Prerequisites

External tools (must be on `$PATH`, same as the original project):

- **R** ≥ 4.0 (only base R — no R packages required)
- **SAMtools** ≥ 1.15
- **BWA** ≥ 0.7.17
- **SPAdes** ≥ 3.13 (only needed when `--flanking` is used)
- **Python 3** ≥ 3.8

Python packages needed by the inner helper scripts (install via pip or the
provided conda yaml):

- `pysam`, `pyfastx`, `numpy`, `pandas`, `matplotlib`

**Also required:**

- A local clone of the **original** `gAIRR-wgs` repository, because its
  `gAIRR_suite/scripts/` and `gAIRR_suite/material/` directories are reused.

## Installation

### Option A — conda (self-contained, recommended on laptops / workstations)

```bash
# 1. Clone both repos into the same parent directory
cd ~/work
git clone https://github.com/AnthonyDNA/gAIRR-wgs.git gAIRR-wgs-main
git clone https://github.com/AnthonyDNA/gAIRR-wgs-R-version.git gAIRR-wgs-R-version

# 2. Create the conda env (same one the original uses; adds r-base)
cd gAIRR-wgs-R-version
conda env create -f gAIRR-wgs_env.yaml
conda activate gAIRR-wgs

# 3. Patch: symlink scripts/ and material/ into the R tree, and install
#    gAIRR_wgs / gAIRR_extract onto $PATH as R shims.
bash clone_script.sh
```

If your original repo lives somewhere unusual, point at it explicitly:

```bash
GAIRR_ORIG=/path/to/gAIRR-wgs bash clone_script.sh
```

### Option B — environment modules (typical on HPC)

If your cluster uses environment modules rather than conda, you can load the
tools from modules and install the Python helper packages into a venv. See
`hpc/taiwania3/` for a fully worked example (Taiwania3-specific, but the
pattern is portable).

## Usage

Exactly like the original. After `clone_script.sh`:

```bash
gAIRR_wgs \
    -wd <work_dir>  \
    -id <sample_ID> \
    -rd1 <read.R1.fastq.gz> \
    -rd2 <read.R2.fastq.gz>
    # -lc defaults to all six loci: TRV TRJ TRD IGV IGJ IGD
```

Or invoke the R file directly (useful if you don't want to install shims):

```bash
Rscript gAIRR_suite_R/gAIRR_wgs.R \
    -wd  target_call \
    -id  HG002 \
    -rd1 HG002.R1.fastq.gz \
    -rd2 HG002.R2.fastq.gz \
    -lc  TRV TRJ TRD IGV IGJ IGD
```

All the flags from the Python CLI are preserved, with identical semantics:

```
-wd/--work_dir   output directory                    [default: target_call/]
-lc/--locus      target locus list                   [default: TRV TRJ TRD IGV IGJ IGD]
-id/--sample_id  sample identifier
-rd1/--read1     pair-end 1 FASTQ                    (required)
-rd2/--read2     pair-end 2 FASTQ                    (required)
    --flanking   also run flanking-sequence analysis
    --keep       keep intermediate .sam files
-t/--thread      thread count                         [default: max]
-REF/--reference path to material/ directory          [default: next to the R script]
```

### gAIRR_extract

```bash
gAIRR_extract \
    -l  <samples.tsv> \
    -o  <output_dir>  \
    [-b <bed>] [-r <ref.fa>] [-t <threads>]
```

`samples.tsv` has the same two supported shapes as the original:

```tsv
# BAM / CRAM
sample_ID   path1
HG001       /data/HG001.cram
HG002       /data/HG002.bam
```

```tsv
# FASTQ
sample_ID   path1                     path2
Sample1     /data/Sample1_R1.fastq.gz /data/Sample1_R2.fastq.gz
```

## Test

The test FASTQs used for validation are the HG002 paired-end reads included with the original repo at `gAIRR-wgs-main/test_sample/HG002_part_R{1,2}.fastq.gz`.

To run the test:

```bash
# 1. copy the FASTQs into this repo's test_sample/ (they aren't committed here
#    to avoid duplicating binary data):
cp ../gAIRR-wgs-main/test_sample/HG002_part_R{1,2}.fastq.gz test_sample/

# 2. edit the four path variables in test_sample/test_script.sh, then run:
cd test_sample
bash test_script.sh
```

For an HPC-friendly test driver (SLURM + modules), use `hpc/taiwania3/01_smoke_test_modules.sbatch` instead.

## What exactly is "R" about this version?

Only the two user-facing orchestrators have been ported:

| Original (Python / shell)                | R version                               |
| ---------------------------------------- | --------------------------------------- |
| `gAIRR_suite/gAIRR_wgs.py`               | `gAIRR_suite_R/gAIRR_wgs.R`             |
| `gAIRR_suite/gAIRR_extract.sh`           | `gAIRR_suite_R/gAIRR_extract.R`         |

The ~50 inner analysis scripts (`parse_cluster_realign.py`, `novel_allele.sh`, `flanking_sequence.sh`, etc.) are **intentionally reused as-is** from the original repository. Reimplementing them in R would risk introducing silent numerical or parsing differences from the published pipeline, which is exactly what we wanted to avoid. The R entry points call into those scripts through `system()`, just as the Python version did via `subprocess.call(..., shell=True)`.

## Parity with the Python version

The R port is designed to be bit-for-bit equivalent to the Python version
because it delegates all sequence work to the same underlying scripts.

**Reference validation on NCHC Taiwania3** using the HG002_part test FASTQs across all six gAIRR loci (TCRV, TCRJ, TCRD_plusHep, BCRV, BCRJ, BCRD_plusHep):

- 87 output files in identical directory structures (both versions produce the same set of files, nothing extra, nothing missing).
- All 36 deterministic text outputs match **md5-identically** between the Python and R pipelines: final allele-call reports (`gAIRR-call_report.rpt`, `read_depth_calling_by_bwa.rpt`), novel-allele FASTAs (`*_with_novel.fasta`), corrected-allele chains (`corrected_alleles_{raw,filtered,extended}.fasta`), and BWA coverage
  reports (`bwa_alleles_all.txt`, `bwa_alleles_cover_link.rpt`).

To reproduce the parity check yourself, see `hpc/taiwania3/04_parity_check_modules.sbatch`.

## Notes

- **Output compatibility.** Every output file is produced by the same underlying tools (bwa, samtools, SPAdes, plus the `scripts/*.py` helpers), so the R version is byte-identical to the original's output for identical inputs. See the parity-check section above for evidence.
- **`system()` vs `system2()`.** The R wrappers deliberately use `system()` (shell invocation) rather than `system2()` because a handful of the original cleanup commands use glob expansion (e.g., `rm .../*.sam`) that only a shell can perform. This exactly matches the Python version's `subprocess.call(..., shell=True)`.
- **No R package dependencies.** Argument parsing, OS detection, and logging are all written in base R so `Rscript` + a working `R` install is enough.