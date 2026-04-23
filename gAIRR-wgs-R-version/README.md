# gAIRR-wgs (R version)

An R port of the user-facing entry points of [gAIRR-wgs](https://github.com/AnthonyDNA/gAIRR-wgs). The original project is a bioinformatics pipeline that genotypes germline **T-cell receptor (gTR)** genes from standard whole-genome sequencing data (paired-end 150 bp, ≥30× coverage).

## Prerequisites

External tools (must be on `$PATH`, same as the original project):

- **R** ≥ 4.0
- **SAMtools** ≥ 1.15
- **BWA** ≥ 0.7.17
- **SPAdes** ≥ 3.13
- **Python 3** ≥ 3.7.13

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
git clone https://github.com/AnthonyDNA/gAIRR-wgs.git gAIRR-wgs-main

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
    # -lc defaults to all six loci: TRV TRJ TRD
```

Or invoke the R file directly (useful if you don't want to install shims):

```bash
Rscript gAIRR_suite_R/gAIRR_wgs.R \
    -wd  target_call \
    -id  HG002 \
    -rd1 HG002.R1.fastq.gz \
    -rd2 HG002.R2.fastq.gz \
    -lc  TRV TRJ TRD
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

The test FASTQs used for validation are the HG002 paired-end reads included with the original repo at `gAIRR-wgs/test_sample/HG002_part_R{1,2}.fastq.gz`.

To run the test:

```bash
# 1. copy the FASTQs into this repo's test_sample/ (they aren't committed here
#    to avoid duplicating binary data):
cp ../gAIRR-wgs/test_sample/HG002_part_R{1,2}.fastq.gz test_sample/

# 2. edit the four path variables in test_sample/test_script.sh, then run:
cd test_sample
bash test_script.sh
```

For an HPC-friendly test driver (SLURM + modules), use `hpc/taiwania3/01_smoke_test_modules.sbatch` instead.
