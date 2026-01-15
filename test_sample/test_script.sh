word_dir=path/to/your/work_dir
sample_ID=HG002
read1=path/to/your/HG002.R1.fastq.gz
read2=path/to/your/HG002.R2.fastq.gz

gAIRR_wgs \
    -wd $word_dir \
    -id $sample_ID \
    -rd1 $read1 \
    -rd2 $read2 \
    -lc TRV TRJ TRD