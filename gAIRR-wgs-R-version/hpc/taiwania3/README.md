# Example deployment: NCHC Taiwania3 (modules)

This folder contains a concrete, reproducible example of running the R port of gAIRR-wgs on [NCHC Taiwania3](https://www.nchc.org.tw/) using the cluster's **environment-module** system (no conda). It is *one* example deployment — if you're on a different cluster the same structure should work with minor
module-name edits. All scripts use the SBATCH header shape commonly used on Taiwania3:

```
#SBATCH -A CHANGE_ME_SLURM_ACCOUNT     
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
```

Before submitting any of the scripts, edit the `-A` line to your project
account. Email notifications are commented out by default; uncomment the two
`--mail-user` / `--mail-type` lines at the top of each script if you want
them.

## Adapting to a different HPC

The only Taiwania3-specific bits are:

- the `module load ...` block (replace with your cluster's module names),
- the `-p ngs92G` partition and the `-A ...` account,
- the `/staging/biology/$USER` path (replace with your cluster's scratch
  / project directory).

Everything downstream — the R wrapper invocation, the venv activation, the
inner pipeline — is portable.

## Install layout (the scripts assume this tree)

```
/staging/biology/$USER/
├── gAIRR-wgs-R-version/       # this R port (you put it here yourself)
├── gAIRR-wgs-main/            # original Python/shell repo (00_setup_modules.sh clones it for you)
├── gAIRR-venv/                # python venv with pysam / pyfastx / numpy / pandas / matplotlib
└── reference/                 # (optional) hg38_mod.fa for gAIRR_extract
```

## Phase 0 — one-time setup (runs on the login node)

```bash
cd /staging/biology/$USER/gAIRR-wgs-R-version/hpc/taiwania3
bash 00_setup_modules.sh
```

What it does:

1. `module load biology`, `R/4.4.1`, `BWA/0.7.17`, `Samtools/1.15.1`, `python/3.12.2`, `SPAdes/3.15.5`, then prints the resolved tool paths.
2. `git clone` the original `gAIRR-wgs` repo next to the R version (if missing) — we need its `scripts/` and `material/` folders.
3. Create a Python venv at `/staging/biology/$USER/gAIRR-venv` on top of `python/3.12.2`.
4. `pip install` the helper-script deps (pysam, pyfastx, numpy, pandas, matplotlib).
5. Run `clone_script.sh` to symlink `scripts/` and `material/` into the R tree.
6. Print a summary sanity check: `Rscript gAIRR_wgs.R --help` and a Python `import` test.

## Phase 1 — smoke test on the included test FASTQs 

Stage the 300 KB test FASTQs that ship with the original repo (or your own
small pair) at:

```
/staging/biology/$USER/gAIRR-wgs-R-version/test_sample/HG002_part_R1.fastq.gz
/staging/biology/$USER/gAIRR-wgs-R-version/test_sample/HG002_part_R2.fastq.gz
```

Then submit:

```bash
cd /staging/biology/$USER/gAIRR-wgs-R-version/hpc/taiwania3
sbatch 01_smoke_test_modules.sbatch

squeue -u $USER        # watch it
tail -f out_smoke.log  # follow progress
```

The job loads the six modules, activates the venv, and runs the R entry script directly. Pass = the R port dispatches each requested locus end-to-end without error. Expected outputs (for TRV; other loci produce analogous folders):

- `smoke_out/HG002_smoke/HG002_smoke_TCRV/` — alignment + allele-call
  intermediates
- `smoke_out/HG002_smoke/HG002_smoke_TCRV_novel/TCRV_with_novel.fasta` —
  novel-allele FASTA
- job exits 0 with `=== Smoke test PASSED ===`

## Phase 2 — gAIRR_wgs on your real WGS FASTQs

Edit the paths at the top of `02_run_gAIRR_wgs_modules.sbatch`, then:

```bash
sbatch 02_run_gAIRR_wgs_modules.sbatch
```

Defaults to all six loci (TRV TRJ TRD IGV IGJ IGD). Change `LOCI=...` for a subset. Set `FLANKING="--flanking"` for flanking analysis (needs SPAdes; longer runtime). 

## Phase 3 — gAIRR_extract on BAM / FASTQ samples

`gAIRR_extract` needs:

- a TSV sample list (see `sample_list_template.tsv`),
- a BED file that ships with the original repo at `gAIRR-wgs-main/gAIRR_suite/material/extract_materials/gAIRR_allele.bed`,
- an **hg38 reference FASTA** (not in the repo) with a BWA index alongside.

```bash
# 1. prepare your sample list
cp sample_list_template.tsv samples.tsv
vim samples.tsv     # fill in sample_ID + BAM path (or paired FASTQ paths)

# 2. make sure REF_FA in the script points at an hg38 .fa that already has
#    .bwt / .pac / .sa / .amb / .ann next to it.  If not, run once on a login
#    node with modules loaded:
#       module load biology BWA/0.7.17 Samtools/1.15.1
#       bwa index  /path/to/hg38.fa
#       samtools faidx /path/to/hg38.fa

# 3. submit
sbatch 03_run_gAIRR_extract_modules.sbatch
```

## Phase 4 — parity check (R vs. Python on the same sample)

The point of the port is to produce the same results. `04_parity_check_modules.sbatch` runs the original Python version and the R version on the same input, then diffs final FASTAs and allele-call text reports. A zero-diff result is the evidence that the R wrapper is a faithful
port.

```bash
sbatch 04_parity_check_modules.sbatch
```

The reference parity run on NCHC Taiwania3 with the HG002_part test FASTQs produced **87 output files in identical directory structures** and all **36 deterministic text outputs match md5-identically** across all six gAIRR loci (TCRV, TCRJ, TCRD_plusHep, BCRV, BCRJ, BCRD_plusHep).

## A note on FASTQ vs BAM input

| Tool            | FASTQ | BAM/CRAM |
| --------------- | :---: | :------: |
| `gAIRR_wgs`     |  required |  no  |
| `gAIRR_extract` |  yes  |  yes     |

If all you have is a BAM, run Phase 3 first to produce AIRR-filtered +
unmapped FASTQs (`on_AIRR_unmapped_R1.fastq`, `on_AIRR_unmapped_R2.fastq`),
then feed *those* into Phase 2.