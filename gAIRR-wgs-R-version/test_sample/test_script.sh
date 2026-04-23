#!/bin/bash
# Test driver for the R version of gAIRR-wgs.
# Edit the four variables below before running.

work_dir=path/to/your/work_dir
sample_ID=HG002
read1=path/to/your/HG002.R1.fastq.gz
read2=path/to/your/HG002.R2.fastq.gz

# After running clone_script.sh, `gAIRR_wgs` is an R shim on $PATH.
gAIRR_wgs \
    -wd  "$work_dir"  \
    -id  "$sample_ID" \
    -rd1 "$read1"     \
    -rd2 "$read2"     \
    -lc  TRV TRJ TRD

# If you prefer to skip the shim, invoke the R script directly:
# Rscript ../gAIRR_suite_R/gAIRR_wgs.R \
#     -wd  "$work_dir" \
#     -id  "$sample_ID" \
#     -rd1 "$read1" \
#     -rd2 "$read2" \
#     -lc  TRV TRJ TRD
