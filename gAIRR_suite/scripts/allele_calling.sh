# path parameters
outer_dir=$1/$4_$2/
allele_name=$2
allele_path=$3
read_path_1=$5
read_path_2=$6
thread=$7
coverage_thrsd=80

# get the script directory
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

db_allele_flanking_path=${script_dir}/../material/${allele_name}_alleles_parsed.fasta
db_allele_original_path=${script_dir}/../material/ori_material/${allele_name}_alleles_parsed.fasta



# Debug output
echo "Allele Path: ${allele_path}"
echo "Read Path 1: ${read_path_1}"
echo "Read Path 2: ${read_path_2}"
echo "Thread count: ${thread}"
echo "Outer directory: ${outer_dir}"
# setting for the data
mkdir -p ${outer_dir}
bwa index ${allele_path}
if [ ${allele_name} == "TCRD_plusHep" ] || [ ${allele_name} == "BCRD_plusHep" ]; then
    echo "[gAIRRCall: gAIRR-wgs] Adjust BWA parameters for shorter alleles..."
    bwa mem -t ${thread} -T 10 ${allele_path} ${read_path_1} ${read_path_2} > ${outer_dir}bwa_read_to_allele_all.sam
else
    bwa mem -t ${thread} ${allele_path} ${read_path_1} ${read_path_2} > ${outer_dir}bwa_read_to_allele_all.sam
fi
echo "[gAIRRCall: gAIRR-wgs] [Filtering]"
############################################
###              Sam to Bam              ###
############################################
echo "sam to bam"
samtools view -b -h ${outer_dir}bwa_read_to_allele_all.sam > ${outer_dir}bwa_read_to_allele_all.bam
if [ $? -ne 0 ]; then
    echo "Error: samtools view -b -h s>b failed"
fi
echo "bam: filter out reads without reference"
############################################
###       Filter out unmapped read       ###
############################################
samtools view -F 4 -h -b ${outer_dir}bwa_read_to_allele_all.bam > ${outer_dir}filtered_bwa_read_to_allele_all.bam
if [ $? -ne 0 ]; then
    echo "Error: samtools view -F 4 b>b failed"
fi
rm ${outer_dir}bwa_read_to_allele_all.bam
############################################
###              Sort Bam                ###
############################################
echo "sort bam"
samtools sort ${outer_dir}filtered_bwa_read_to_allele_all.bam -o ${outer_dir}sorted_filtered_bwa_read_to_allele_all.bam
if [ $? -ne 0 ]; then
    echo "Error: samtools sort failed"
fi
samtools index ${outer_dir}sorted_filtered_bwa_read_to_allele_all.bam
rm ${outer_dir}filtered_bwa_read_to_allele_all.bam
############################################
###              Bam to Sam              ###
############################################
echo "bam to sam"
samtools view -h ${outer_dir}sorted_filtered_bwa_read_to_allele_all.bam > ${outer_dir}bwa_read_to_allele_all.sam
if [ $? -ne 0 ]; then
    echo "Error: samtools view -h b>s failed"
fi

######################################################
###  Extract primary aligned pair-end information  ###
######################################################
samtools view -h -f 0x2 -F 0x104 ${outer_dir}sorted_filtered_bwa_read_to_allele_all.bam > ${outer_dir}primary_pair.sam

# start analysis
echo "[gAIRRCall: gAIRR-wgs] [ALLELE CALLING] Calling alleles..."
python3 ${script_dir}/analyze_read_depth_with_bwa.py -fs  ${outer_dir}bwa_read_to_allele_all.sam \
                                                     -an  ${allele_name} \
                                                     -fa  ${allele_path} \
                                                     -t   ${coverage_thrsd} \
                                                     -foc ${outer_dir}read_depth_calling_by_bwa.rpt \
                                                     -fop ${outer_dir}allele_support_reads.pickle \
                                                     -fdbfnk ${db_allele_flanking_path} \
                                                     -fdbraw ${db_allele_original_path} \
                                                     -pap ${outer_dir}primary_pair.sam \
                                                     #-fv  ${annotation_path}
                                                     #-fs  ${outer_dir}bwa_read_to_allele_all.sam \

python3 ${script_dir}/calling_threshold.py -dp ${outer_dir}read_depth_calling_by_bwa.rpt > ${outer_dir}gAIRR-call_report.rpt
echo "[gAIRRCall: gAIRR-wgs] [ALLELE CALLING] Finished!"
