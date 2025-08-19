import pysam
import argparse
from utils import get_reverse_complement


def parse_fasta(fn_fasta):
    dict_fasta = {}
    f = open(fn_fasta)
    while True:
        allele_name = f.readline().strip()[1:]
        sequence    = f.readline().strip().lower()
        if not allele_name:
            break
        allele_name = allele_name.split('|')[1]
        
        if dict_fasta.get(allele_name):
            print("Error Format: Duiplicated Sequence:", allele_name)
        dict_fasta[allele_name] = sequence
    return dict_fasta


def seq_contain(seq, target_seq):
    if seq == target_seq or seq == get_reverse_complement(target_seq):
        return 1
    elif seq in target_seq or seq in get_reverse_complement(target_seq):
        return 2
    elif target_seq in seq or get_reverse_complement(target_seq) in seq:
        return 3
    else:
        return 0




def self_compare(dict_fasta):
    for allele_name, sequence in sorted(dict_fasta.items()):
        for target_name, target_seq in sorted(dict_fasta.items()):
            if allele_name != target_name:
                contain_flag = seq_contain(sequence, target_seq)
                if contain_flag == 0:
                    pass
                elif contain_flag == 1:
                    print(allele_name, 'is equal to', target_name)
                elif contain_flag == 2:
                    print(allele_name, 'is contain in', target_name)
                else:
                    print(allele_name, 'contains', target_name)





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-fa', '--fasta_file', help='the fasta file checking for duplications')
    args = parser.parse_args()

    fn_fasta = args.fasta_file
    
    dict_fasta = parse_fasta(fn_fasta)
    self_compare(dict_fasta)

