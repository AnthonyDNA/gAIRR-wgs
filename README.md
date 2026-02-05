# gAIRR-wgs
**gAIRR-wgs** is designed for genotyping germline ***T cell receptor*** genes (g*TR*)  from standard whole-genome sequencing (WGS) data (paired-end 150 bp, ~30× coverage or higher).
For targeted panel sequencing, please refer to [gAIRR-call](https://github.com/maojanlin/gAIRRsuite).

## Prerequisite programs:
* SAMTOOLS (1.18)
* BWA aligner (0.7.17)
* SPAdes assembler (v3.13.0)

## Installation (recommend)
1. Create work directory for gAIRR-wgs
```
mkdir {your_work_dir}
```
```
cd {your_work_dir}
```
2. Clone gAIRR-wgs script

```
git clone git@github.com:AnthonyDNA/gAIRR-wgs.git
```
3. Create environment management system
```
conda env create -f gAIRR-wgs_env.yaml
```
4. Activate virtual environment
```
conda activate gAIRR-wgs
```
5. Patch the gAIRR-wgs module
```
bash clone_script.sh
```

## Usage
```
$ gAIRR_wgs \
    -wd <work_dir> \
    -id <sample_ID> \
    -rd1 <read.R1.fastq.gz> \
    -rd2 <read.R2.fastq.gz> \
    -lc <TRV TRJ TRD>
```

## Example
This repository provides test materials under `test_sample/`.

### Test inputs
- `test_sample/HG002_part_R1.fastq.gz`
- `test_sample/HG002_part_R2.fastq.gz`

### Run the test
1. Open and modify `test_sample/test_script.sh` as needed.
2. Execute:
   ```
   cd test_sample
   bash test_script.sh
   ```


## gAIRR_extract - Read Extraction Tool
Extract gAIRR reads from BAM/CRAM/FASTQ files

### Usage

```
$ gAIRR_extract \
    -l <samples.tsv> \
    -o <output_dir> \
    [-t <thread>]
```
### Example
#### Sample List Format
**For BAM/CRAM files:**
```tsv
sample_ID	path1
HG001	/path/to/HG001.cram
HG002	/path/to/HG002.bam
```

**For FASTQ files:**
```tsv
sample_ID	path1	path2
Sample1	/path/to/Sample1_R1.fastq.gz	/path/to/Sample1_R2.fastq.gz
Sample2	/path/to/Sample2_R1.fq	/path/to/Sample2_R2.fq
```
#### Command
Execute:
   ```
   gAIRR_extract \
    -l /path/to/sample_list.tsv
    -o /path/to/out_put_dir
    -t 20
   ```
