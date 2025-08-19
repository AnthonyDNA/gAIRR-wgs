import argparse
import pickle
import os
import numpy as np
import sys

# make sure the package modules is in the path
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.dirname(__file__)+'/..')

from utils import get_reverse_complement
from filter_corrected_alleles import parse_perfect_sam, parse_fasta
from parse_contig_realign import parse_CIGAR
import sys

def clean_dm_suffix(seq_name):
    """Remove _DM_ suffix from sequence name if present"""
    if seq_name.endswith('_DM_'):
        return seq_name[:-4]  # Remove the last 4 characters '_DM_'
    return seq_name

def process_novel_allele_name(name):
    """Process novel allele name based on whether it contains _DM_ suffix"""
    if '_DM_' in name:
        # Remove _DM_ suffix and format as >new|TRAV25*01_JW7E|Homo
        cleaned_name = clean_dm_suffix(name)
        return cleaned_name
    else:
        # Keep original format but wrap with |name|
        return "|" + name + "|"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fa', '--fn_alleles',
        help = 'input original allele fasta file'
    )
    parser.add_argument(
        '-fn', '--fn_novel',
        help = 'input fasta file with novel alleles'
    )
    
    parser.add_argument(
        '-fom', '--fo_merged_fasta',
        help = 'output merged fasta file'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



if __name__ == '__main__':
    args = parse_args()
    fn_alleles = args.fn_alleles
    fn_novel   = args.fn_novel
    fo_merged_fasta = args.fo_merged_fasta

    dict_original_allele = parse_fasta(fn_alleles)
    dict_novel_allele    = parse_fasta(fn_novel)

    f_om = open(fo_merged_fasta, 'w')
    for name in sorted(dict_original_allele.keys()):
        f_om.write(">" + name.split()[0] + "\n")
        f_om.write(dict_original_allele[name] + "\n")
    if dict_novel_allele == {'':''}: #empty novel results
        pass
    else:
        for name in sorted(dict_novel_allele.keys()):
            if ("extend" in name) == False:
                processed_name = process_novel_allele_name(name)
                f_om.write(">" + processed_name + "\n")
                f_om.write(dict_novel_allele[name] + "\n")
    f_om.close()

