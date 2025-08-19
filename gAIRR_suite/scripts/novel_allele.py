import subprocess
import sys
import os
import argparse
import io
from contextlib import redirect_stdout

import parse_cluster_realign


def main(arguments=None):
    # ${workspace}/${person_name} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
    parser = argparse.ArgumentParser(description="gAIRR-call pipeline of gAIRR-suite.")
    parser.add_argument('-wd', '--work_dir', help="Path to output directory ['target_call/'].", default="target_call/")
    parser.add_argument('-lc', '--locus', help="Target Locus [TRV TRJ TRD IGV IGJ IGD]", nargs='+', default=['TRV', 'TRJ', 'TRD', 'IGV', 'IGJ', 'IGD'])
    parser.add_argument('-id', '--sample_id', help="Sample ID ['sample']")
    parser.add_argument('-rd1', '--read1', help="Path to gAIRR-seq. Pair-end 1", required=True)
    parser.add_argument('-rd2', '--read2', help="Path to gAIRR-seq. Pair-end 2", required=True)
    parser.add_argument('--flanking', help="Option to do flanking sequence analysis.", action='store_true')
    parser.add_argument('-t', '--thread', help="Number of threads to use [max].", type=int)
    args = parser.parse_args(arguments)

    workspace     = args.work_dir
    target_locus  = args.locus
    person_name   = args.sample_id
    path_read1    = args.read1
    path_read2    = args.read2
    flag_flanking = args.flanking
    thread        = args.thread
    if thread == None:
        thread = get_max_thread()
    path_module   = os.path.dirname(__file__) + '/'
    path_material = path_module + '../example/material/'
    
    



    # Index and align the gAIRR-seqs
    out_dir = workspace+'/'+person_name+'_'+allele_name+'_novel/'
    subprocess.call("mkdir -p " + workspace, shell=True)
    subprocess.call("bwa", "index", path_allele)
    if 'plusHep' in allele_name: # D genes, apply special bwa parameters
        command_prefix = "bwa mem -a -T 10 -t"
    else:
        command_prefix = "bwa mem -a -t"
    command = (command_prefix, str(thread), path_allele, read_1, read_2, '>', outer_dir+"/bwa_read_to_allele.sam")
    subprocess.call(' '.join(command), shell=True)

    # Analyze the reads
    print("[gAIRR-call] [NOVEL ALLELE] Finding novel alleles...")
    if os.path.isfile(out_dir + 'corrected_alleles_raw.fasta') == True:
        subprocess.call("rm " + out_dir + 'corrected_alleles_raw.fasta', shell=True)
    command =  ['-fs',  outer_dir+'bwa_read_to_allele.sam', \
                '-fc',  allele_path, \
                '-for', 'test.txt', \
                '-foc', outer_dir+'corrected_alleles_raw.fasta'] \
                #>    ${outer_dir}bwa_alleles_cover_link.rpt \
                #2>   ${outer_dir}bwa_alleles_all.txt
    
    f = io.StringIO()
    with redirect_stdout(f):
        parse_cluster_realign.main(command)
    out = f.getvalue()
    f = open(outer_dir+'bwa_alleles_cover_link.rpt', 'w')
    f.write(out)










    dict_locus_map = {'TRV':'TCRV', 'TRJ':'TCRJ', 'TRD':'TCRD_plusHep', 'IGV':'BCRV', 'IGJ':'BCRJ', 'IGD':'BCRD_plusHep'}
    for locus_name in target_locus:
        if dict_locus_map.get(locus_name) == None:
            print('Only TRV, TRJ, TRD, IGV, IGJ, IGD are allowed for target locus call.')
            parser.print_usage()
            exit(2)
    list_allele_names = [dict_locus_map[locus_name] for locus_name in target_locus]
    

    subprocess.call("mkdir -p " + workspace, shell=True)
    subprocess.call("mkdir -p " + workspace+'/'+person_name, shell=True)
    # call novel alleles
    for allele_name in list_allele_names:
        allele_path = path_material+allele_name+"_alleles_parsed.fasta"
        print("[gAIRR-call]", person_name, allele_name, "novel allele calling...")
        #./scripts/novel_allele.sh ${workspace}/${person_name} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
        command = ' '.join(["bash", path_module+"novel_allele.sh", workspace+'/'+person_name, allele_name, allele_path, person_name, path_read_1, path_read_2])
        subprocess.call(command, shell=True)
        
        allele_path = workspace+'/'+person_name+'/'+person_name+'_'+allele_name+"_novel/"+allele_name+"_with_novel.fasta"
        print("[gAIRR-call]", person_name, allele_name, "allele calling...")
        #./scripts/allele_calling.sh ${workspace}/${person_name} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2}
        command = ' '.join(["bash", path_module+"allele_calling.sh", workspace+'/'+person_name, allele_name, allele_path, person_name, path_read_1, path_read_2])
        subprocess.call(command, shell=True)
        
        if flag_flanking:
            print("[gAIRR-call]", person_name, allele_name, "flanking sequence calling...")
            #./scripts/flanking_sequence.sh ${workspace}/${person_name} ${allele_name} ${allele_path} ${person_name} ${read_path_1} ${read_path_2} ${path_SPAdes}
            command = ' '.join(["bash", path_module+"allele_calling.sh", workspace+'/'+person_name, allele_name, allele_path, person_name, path_read_1, path_read_2, 'spades.py'])
            subprocess.call(command, shell=True)

    




if __name__ == "__main__":
    main()


