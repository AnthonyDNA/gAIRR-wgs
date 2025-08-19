import argparse
import pickle
import os
import numpy as np

import sys

# make sure the package modules is in the path
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.dirname(__file__)+'/..')

from scripts.utils import get_reverse_complement
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fa', '--fn_original_alleles',
        help = 'input original allele fasta file'
    )
    parser.add_argument(
        '-fca', '--fn_corrected_alleles',
        help = 'input corrected allele fasta file'
    )
    
    parser.add_argument(
        '-fof', '--fo_filtered_alleles',
        help = 'output filtered corrected allele file'
    )
    parser.add_argument(
        '-foe', '--fo_extended_alleles',
        help = 'output extended allele file'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_fasta(fn_fasta):
    '''parse the fasta file into a dictionary'''
    # dict_name_SEQ {}
    #  - keys: seq_name
    #  - values: seq_SEQ
    dict_name_SEQ = {}
    with open(fn_fasta, 'r') as f_f:
        seq_name = ""
        seq_SEQ = ""
        for line in f_f:
            if line[0] == '>':
                if seq_name != "":
                    if dict_name_SEQ.get(seq_name):
                        print("WARNING! Duplicate sequence name:", seq_name)
                        new_name = assign_new_name_basic(seq_name, dict_name_SEQ)
                        dict_name_SEQ[new_name] = seq_SEQ
                    else:
                        dict_name_SEQ[seq_name] = seq_SEQ
                seq_name = line.strip()[1:]
                try:
                    seq_name = seq_name.split()[0]
                except:
                    pass
                seq_SEQ = ""
            else:
                seq_SEQ += line.strip()
        if dict_name_SEQ.get(seq_name):
            new_name = assign_new_name_basic(seq_name, dict_name_SEQ)
            dict_name_SEQ[new_name] = seq_SEQ
            print("WARNING! Duplicate sequence name:", seq_name)
        else:
            dict_name_SEQ[seq_name] = seq_SEQ
    return dict_name_SEQ


def parse_perfect_sam(fn_sam):
    list_perfect_fields = []
    list_mismatch_fields = []
    with open(fn_sam, 'r') as f_s:
        for line in f_s:
            if line[0] != '@':
                fields = line.strip().split()
                if "NM:i:0" in line and ('S' in fields[5] or 'H' in fields[5]) == False:
                    list_perfect_fields.append(fields)
                elif fields[2] != '*':
                    list_mismatch_fields.append(fields)
    return list_perfect_fields, list_mismatch_fields


def assign_new_name_basic(basic_name, dict_target):
    ext_num = 0
    while True:
        tmp_name = basic_name + '_' + str(ext_num)
        if dict_target.get(tmp_name):
            ext_num+=1
        else:
            return tmp_name


def assign_new_name(basic_name, str_interval, dict_target):
    if '|' in basic_name:
        try:
            basic_name = basic_name.split('|')[1]
        except:
            pass
    ext_num = 0
    while True:
        tmp_name = basic_name + str_interval + str(ext_num)
        if dict_target.get(tmp_name):
            ext_num += 1
        else:
            return tmp_name


def assign_new_name_preserve_dm(basic_name, str_interval, dict_target):
    """For _DM_ sequences, preserve the full name structure without extracting middle part"""
    ext_num = 0
    while True:
        tmp_name = basic_name + str_interval + str(ext_num)
        if dict_target.get(tmp_name):
            ext_num += 1
        else:
            return tmp_name


def clean_dm_suffix(seq_name):
    """Remove _DM_ suffix from sequence name if present"""
    if seq_name.endswith('_DM_'):
        return seq_name[:-4]  # Remove the last 4 characters '_DM_'
    return seq_name


def duplicate_trim_set_with_2nd_set(dict_target_allele_SEQ, dict_fixed_allele_SEQ, ext_flag=True, ext_thrd=0.70, ori_flag=False, dm_names_set=None):
    # if any target_SEQ is subseq of any fixed_SEQ, the target_SEQ is popped
    # if ext_flag is True, when a fixed_SEQ is a subseq of a target_SEQ, target_SEQ is kept and target_name is changed into fixed_name_ext
    # the ext_flag is effective if only len(fixed_SEQ) > len(target_SEQ)*ext_thrd
    # dm_names_set: set of cleaned names that were originally _DM_ sequences
    dict_trimmed_allele_SEQ = {}
    if dm_names_set is None:
        dm_names_set = set()
    
    for (t_name, t_SEQ) in dict_target_allele_SEQ.items():
        assign_name = t_name
        
        # Check if this is a _DM_ sequence (either currently has _DM_ suffix or was originally a _DM_ sequence)
        is_dm_sequence = t_name.endswith('_DM_') or t_name in dm_names_set
        
        for (f_name, f_SEQ) in dict_fixed_allele_SEQ.items():
            if t_SEQ.upper() in f_SEQ.upper() or t_SEQ.upper() in get_reverse_complement(f_SEQ.upper()): # t_SEQ is the subseq
                assign_name = False
                break
            elif f_SEQ.upper() in t_SEQ.upper() or f_SEQ.upper() in get_reverse_complement(t_SEQ.upper()): # an f_SEQ is the subseq
                if ext_flag == True and len(f_SEQ) >= len(t_SEQ)*ext_thrd:
                    ext_num = 0
                    if is_dm_sequence:
                        # For _DM_ sequences, use the clean name for extension
                        if t_name.endswith('_DM_'):
                            # clean_name = clean_dm_suffix(t_name)
                            clean_name = t_name
                        else:
                            clean_name = t_name
                        assign_name = assign_new_name_preserve_dm(clean_name, '/extend-', dict_trimmed_allele_SEQ)
                    else:
                        assign_name = assign_new_name(f_name, '/extend-', dict_trimmed_allele_SEQ)
                else:
                    assign_name = False
        
        if assign_name:
            # When ori_flag is True, keep original names
            if ori_flag == True:
                dict_trimmed_allele_SEQ[t_name] = t_SEQ
            else:
                if 'extend' in assign_name:
                    dict_trimmed_allele_SEQ[assign_name] = t_SEQ
                else:
                    # Special handling for _DM_ sequences
                    if is_dm_sequence:
                        # Remove _DM_ suffix and use the clean name directly
                        if t_name.endswith('_DM_'):
                            # clean_name = clean_dm_suffix(assign_name)
                            clean_name = assign_name
                        else:
                            clean_name = assign_name
                        dict_trimmed_allele_SEQ[clean_name] = t_SEQ
                    else:
                        # Regular processing: add /novel- prefix
                        assign_name = assign_new_name(assign_name, '/novel-', dict_trimmed_allele_SEQ)
                        dict_trimmed_allele_SEQ[assign_name] = t_SEQ
    return dict_trimmed_allele_SEQ



if __name__ == "__main__":
    args = parse_args()
    fn_original_alleles = args.fn_original_alleles
    fn_corrected_alleles = args.fn_corrected_alleles
    fo_filtered_alleles = args.fo_filtered_alleles
    fo_extended_alleles = args.fo_extended_alleles

    dict_o_allele_SEQ = parse_fasta(fn_original_alleles)
    dict_c_allele_SEQ = parse_fasta(fn_corrected_alleles)
    
    # Pre-process: Clean all _DM_ suffixes from corrected alleles
    dict_c_allele_SEQ_cleaned = {}
    dm_names_set = set()
    for name, seq in dict_c_allele_SEQ.items():
        if name.endswith('_DM_'):
            # clean_name = clean_dm_suffix(name)
            clean_name = name
            dict_c_allele_SEQ_cleaned[clean_name] = seq
            dm_names_set.add(clean_name)
        else:
            dict_c_allele_SEQ_cleaned[name] = seq
    
    # Use the cleaned dictionary for all subsequent processing
    dict_ref_trimmed_allele_SEQ = duplicate_trim_set_with_2nd_set(dict_c_allele_SEQ_cleaned, dict_o_allele_SEQ, dm_names_set=dm_names_set)
    dict_shrink = {}
    for name, SEQ in dict_ref_trimmed_allele_SEQ.items():
        dict_shrink[name+'_prefix'] = SEQ[:-1]
        dict_shrink[name+'_suffix'] = SEQ[1:]
    dict_self_trimmed_allele_SEQ = duplicate_trim_set_with_2nd_set(dict_ref_trimmed_allele_SEQ, dict_shrink, ext_flag=True, ext_thrd=0, ori_flag=True, dm_names_set=dm_names_set)
    set_SEQ = set()
    for name, SEQ in sorted(dict_self_trimmed_allele_SEQ.items()):
        if SEQ.upper() in set_SEQ:
            dict_self_trimmed_allele_SEQ.pop(name)
        else:
            set_SEQ.add(SEQ.upper())
            set_SEQ.add(get_reverse_complement(SEQ.upper()))

    f_of = open(fo_filtered_alleles, 'w')
    f_oe = open(fo_extended_alleles, 'w')
    for allele_name in sorted(dict_self_trimmed_allele_SEQ.keys()):
        if 'extend' in allele_name:
            f_oe.write(">" + allele_name + '\n')
            f_oe.write(dict_self_trimmed_allele_SEQ[allele_name] + '\n')
        else:
            f_of.write(">" + allele_name + '\n')
            f_of.write(dict_self_trimmed_allele_SEQ[allele_name] + '\n')
    f_of.close()
    f_oe.close()