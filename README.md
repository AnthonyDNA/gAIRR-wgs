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
$ gAIRR_call -wd <work_dir> -id <sample_ID> -rd1 <read.R1.fastq.gz> -rd2 <read.R2.fastq.gz> -lc <TRV TRJ TRD IGV IGJ IGD>
```
