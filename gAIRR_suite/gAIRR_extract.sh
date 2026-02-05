#!/bin/bash
# gAIRR Read Extraction Tool

# Default
DEFAULT_BED="./material/extract_materials/gAIRR_allele.bed"
DEFAULT_REF="./material/references/hg38_mod.fa"
DEFAULT_THREADS=20

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required arguments:
  -l, --list FILE          Sample list TSV file with columns: sample_ID, path1, [path2]
  -o, --outdir DIR         Output directory

Optional arguments:
  -b, --bed FILE           BED file defining AIRR loci regions
                           (default: $DEFAULT_BED)
  -r, --ref FILE           Reference genome FASTA file
                           (default: $DEFAULT_REF)
  -t, --threads INT        Number of threads (default: $DEFAULT_THREADS)
  -h, --help              Show this help message

Sample list format:
  - For BAM/CRAM: sample_ID<TAB>path_to_bam_or_cram
  - For FASTQ:    sample_ID<TAB>R1.fastq<TAB>R2.fastq
  
  The tool auto-detects input format based on file extension.

Example:
  # Using default bed and reference
  $0 -l samples.tsv -o output_dir

  # Using custom bed and reference
  $0 -l samples.tsv -b custom.bed -r custom.fa -o output_dir -t 30

EOF
}

BED_FILE="$DEFAULT_BED"
REF_FASTA="$DEFAULT_REF"
THREADS="$DEFAULT_THREADS"

while [[ $# -gt 0 ]]; do
    case $1 in
        -l|--list)
            SAMPLE_LIST="$2"
            shift 2
            ;;
        -b|--bed)
            BED_FILE="$2"
            shift 2
            ;;
        -r|--ref)
            REF_FASTA="$2"
            shift 2
            ;;
        -o|--outdir)
            BASE_OUTDIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

if [[ -z "$SAMPLE_LIST" || -z "$BASE_OUTDIR" ]]; then
    echo "Error: Missing required arguments (-l and -o are required)"
    print_usage
    exit 1
fi

for file in "$SAMPLE_LIST" "$BED_FILE" "$REF_FASTA"; do
    if [[ ! -f "$file" ]]; then
        echo "Error: File not found: $file"
        exit 1
    fi
done

mkdir -p "$BASE_OUTDIR"

DATE=$(date +%Y%m%d%H%M%S)
LOGFILE="${BASE_OUTDIR}/gAIRR_extraction_${DATE}.log"
exec > >(tee -a "$LOGFILE") 2>&1

# Detect file format
detect_format() {
    local filepath="$1"
    case "$filepath" in
        *.bam) echo "bam" ;;
        *.cram) echo "cram" ;;
        *.fastq.gz|*.fq.gz) echo "fastq.gz" ;;
        *.fastq|*.fq) echo "fastq" ;;
        *) echo "unknown" ;;
    esac
}

# Process BAM/CRAM files
process_alignment() {
    local sample_id="$1"
    local input_path="$2"
    local output_dir="$3"
    
    echo "[${sample_id}] Processing alignment file: $input_path"
    
    # Extract gAIRR reads
    echo "[${sample_id}] Extracting reads on AIRR loci..."
    local extract_start=$(date +%s)
    samtools view -@ "$THREADS" -bh -L "$BED_FILE" -T "$REF_FASTA" "$input_path" | \
        samtools collate -@ "$THREADS" -O -u - | \
        samtools fastq -1 "${output_dir}/on_AIRR_R1.fastq" \
                       -2 "${output_dir}/on_AIRR_R2.fastq" \
                       -0 /dev/null -s /dev/null -n -
    local extract_end=$(date +%s)
    
    # Extract unmapped reads
    echo "[${sample_id}] Extracting unmapped reads..."
    local unmapped_start=$(date +%s)
    samtools view -@ "$THREADS" -bh -f 4 -T "$REF_FASTA" "$input_path" | \
        samtools collate -@ "$THREADS" -O -u - | \
        samtools fastq -1 "${output_dir}/unmapped_R1.fastq" \
                       -2 "${output_dir}/unmapped_R2.fastq" \
                       -0 /dev/null -s /dev/null -n -
    local unmapped_end=$(date +%s)
    
    # Merge
    echo "[${sample_id}] Merging reads..."
    local merge_start=$(date +%s)
    cat "${output_dir}/on_AIRR_R1.fastq" "${output_dir}/unmapped_R1.fastq" > \
        "${output_dir}/on_AIRR_unmapped_R1.fastq"
    cat "${output_dir}/on_AIRR_R2.fastq" "${output_dir}/unmapped_R2.fastq" > \
        "${output_dir}/on_AIRR_unmapped_R2.fastq"
    
    # Cleanup intermediate files
    rm "${output_dir}/on_AIRR_R1.fastq" "${output_dir}/on_AIRR_R2.fastq" \
       "${output_dir}/unmapped_R1.fastq" "${output_dir}/unmapped_R2.fastq"
    local merge_end=$(date +%s)

    echo "  - Read extraction time: $((extract_end - extract_start)) seconds"
    echo "  - Unmapped read extraction time: $((unmapped_end - unmapped_start)) seconds"
    echo "  - Merge time: $((merge_end - merge_start)) seconds"
}

# Process FASTQ files
process_fastq() {
    local sample_id="$1"
    local r1_path="$2"
    local r2_path="$3"
    local output_dir="$4"
    
    echo "[${sample_id}] Processing FASTQ files"
    echo "[${sample_id}] R1: $r1_path"
    echo "[${sample_id}] R2: $r2_path"
    
    # Map reads to reference
    echo "[${sample_id}] Mapping reads to reference..."
    bwa mem -t "$THREADS" -T 10 "$REF_FASTA" "$r1_path" "$r2_path" | \
        samtools view -@ "$THREADS" -bh - | \
        samtools sort -@ "$THREADS" -o "${output_dir}/aligned.bam" -
    
    samtools index -@ "$THREADS" "${output_dir}/aligned.bam"
    
    # Process as alignment file
    process_alignment "$sample_id" "${output_dir}/aligned.bam" "$output_dir"
    
    # Cleanup
    rm "${output_dir}/aligned.bam" "${output_dir}/aligned.bam.bai"
}

# Main processing loop
line_num=0
while IFS=$'\t' read -r sample_id path1 path2; do
    line_num=$((line_num + 1))
    
    # Skip header or empty lines
    if [[ "$sample_id" == "sample_ID" ]] || [[ -z "$sample_id" ]]; then
        continue
    fi
    
    echo ""
    echo "=================================================="
    echo "Processing sample ${line_num}: ${sample_id}"
    echo "=================================================="
    
    # Create sample output directory
    output_dir="${BASE_OUTDIR}/${sample_id}"
    mkdir -p "$output_dir"
    
    start_time=$(date +%s)
    
    # Detect input format
    format=$(detect_format "$path1")
    
    case "$format" in
        bam|cram)
            if [[ ! -f "$path1" ]]; then
                echo "Error: File not found: $path1"
                continue
            fi
            process_alignment "$sample_id" "$path1" "$output_dir"
            ;;
        fastq|fastq.gz)
            if [[ -z "$path2" ]]; then
                echo "Error: R2 path required for FASTQ input"
                continue
            fi
            if [[ ! -f "$path1" ]] || [[ ! -f "$path2" ]]; then
                echo "Error: FASTQ file(s) not found"
                continue
            fi
            process_fastq "$sample_id" "$path1" "$path2" "$output_dir"
            ;;
        *)
            echo "Error: Unknown file format for $path1"
            continue
            ;;
    esac
    
    end_time=$(date +%s)
    total_time=$((end_time - start_time))
    hours=$((total_time / 3600))
    minutes=$((total_time % 3600 / 60))
    seconds=$((total_time % 60))
    
    echo "------------------------------------------------"
    printf "[${sample_id}] Total processing time: %02d:%02d:%02d (hh:mm:ss)\n" \
        $hours $minutes $seconds
    echo "------------------------------------------------"
    
done < "$SAMPLE_LIST"

echo "=================================================="
echo "All samples processed successfully!"
echo "Completed at: $(date)"
echo "=================================================="