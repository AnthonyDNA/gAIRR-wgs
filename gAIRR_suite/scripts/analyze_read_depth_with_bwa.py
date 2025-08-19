import argparse
import pickle
import os
import numpy as np
import sys
from gAIRRwgs_package import parse_sam, get_FL_mark_dict, positions_to_check, FL_region_for_check, parse_fasta3
from collections import defaultdict

# make sure the package modules is in the path
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.dirname(__file__)+'/..')

from utils import get_reverse_complement
from filter_corrected_alleles import parse_perfect_sam, parse_fasta
from parse_contig_realign import parse_CIGAR

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fa', '--fn_alleles',
        help = 'input target allele fasta file'
    )
    parser.add_argument(
        '-fs', '--fn_sam',
        help = 'input short reads to allele sam file'
    )
    parser.add_argument(
        '-an', '--allele_name',
        help = 'allele name'
    )
    parser.add_argument(
        '-t', '--thrsd', type=int,
        help = 'mininum coverage requirement for reads in read depth analysis'
    )
    
    parser.add_argument(
        '-foc', '--fo_calling_report',
        help = 'output calling report file'
    )
    parser.add_argument(
        '-fop', '--fo_grouping_pickle',
        help = 'output allele-read group pickle file'
    )
    parser.add_argument(
        '-fv', '--fn_verify_annotation',
        help = 'input annotation file for verification'
    )
    parser.add_argument(
        '-fdbfnk', '--fn_database_with_flanking', 
        help='input database_with_flanking.fasta'
    )
    parser.add_argument(
        '-fdbraw', '--fn_database_without_flanking', 
        help='input database_without_flanking.fasta'
    )
    parser.add_argument(
        '-pap', '--fn_primary_aligned_pair',
        help='input primary_pair.sam'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parse_pair_sam(primary_sam):
    pair_end_dict = {}
    with open(primary_sam, 'r') as f_s:
        for line in f_s:
            if line[0] != '@':
                fields = line.strip().split()
                read_name = fields[0]
                sum_flag = fields[1]
                ref_name = fields[2]
                read_unique_name = f'{read_name}_{sum_flag}_{ref_name}'
                if read_name not in pair_end_dict:
                    pair_end_dict[read_name]=[]
                pair_end_dict[read_name].append(fields)
    return pair_end_dict

def process_pair(dict_allele_histogram, field1, field2, weight):
    # field1 and field2 are tuples
    allele_name1, allele_name2 = field1[2], field2[2]
    if allele_name1 != allele_name2:
        return True
    else:
        num_list1, operate_list1 = parse_CIGAR(field1[5])
        num_list2, operate_list2 = parse_CIGAR(field2[5])
        start_pos1, start_pos2 = int(field1[3]), int(field2[3])
        end_pos1 = start_pos1 + sum(num for num, op in zip(num_list1, operate_list1) if op == 'M')
        end_pos2 = start_pos2 + sum(num for num, op in zip(num_list2, operate_list2) if op == 'M')
        min_pos = min(start_pos1-1, start_pos2-1)
        max_pos = max(end_pos1-1, end_pos2-1)
        dict_allele_histogram[allele_name1][min_pos:max_pos] += weight
        return False


def parse_perfect_sam_with_S(fn_sam):
    list_perfect_fields = []
    list_mismatch_fields = []

    count_read_cluster_dict = {} 

    with open(fn_sam, 'r') as f_s:
        for line in f_s:
            if line[0] != '@':
                fields = line.strip().split()
                if "NM:i:0" in line:
                    list_perfect_fields.append(fields)
                else:
                    list_mismatch_fields.append(fields)
                
    return list_perfect_fields, list_mismatch_fields

def histogram_read_depth(dict_allele_histogram, list_perfect_fields, dict_allele_reads, salvaged_list, thrsd=100, Pri=False):
    weight_dict = {0: 1.0, 1: 0.5, 2: 0.1, 3: 0, 4: 0, 5: 0} if not Pri else {0: 1.0, 1: 1, 2: 1, 3: 0.5, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0}

    for fields in list_perfect_fields:
        allele_name = fields[2]
        read_name = fields[0]
        start_pos = int(fields[3])
        cigar = fields[5]
        ref_end_pos = len(dict_allele_histogram[allele_name]) + 1
        min_coverage_thrsd = min(thrsd, ref_end_pos-1)
        number, operate = parse_CIGAR(cigar)
        if len(operate) > 3:
            # eprint("Error on CIGAR!", operate)
            continue
        M_idx = None
        for idx, ele in enumerate(operate):
            if ele == 'M':
                M_idx = idx
        if number[M_idx] < min_coverage_thrsd:
            continue
        end_pos = number[M_idx] + start_pos
        if start_pos == 1: # read map to the left side
            if operate[-1] == 'M' or end_pos == ref_end_pos: # right side is open or till the end
                dict_allele_histogram[allele_name][start_pos-1:end_pos-1] += 1
        elif end_pos == ref_end_pos: # read map to the left side
            if operate[0] == 'M': # left side is open
                dict_allele_histogram[allele_name][start_pos-1:end_pos-1] += 1
        elif operate == ['M']: # read is including in the allele
            dict_allele_histogram[allele_name][start_pos-1:end_pos-1] += 1

        dict_allele_reads[allele_name].add(read_name)

    
    for fields, flanking_variant_count in salvaged_list:
        allele_name = fields[2]
        read_name = fields[0]
        start_pos = int(fields[3])
        cigar = fields[5]
        ref_end_pos = len(dict_allele_histogram[allele_name]) + 1
        min_coverage_thrsd = min(thrsd, ref_end_pos-1)

        number, operate = parse_CIGAR(cigar)

        M_idx = None
        for idx, ele in enumerate(operate):
            if ele == 'M':
                M_idx = idx
        if number[M_idx] < min_coverage_thrsd:
            continue
        end_pos = number[M_idx] + start_pos

        weight = weight_dict.get(flanking_variant_count, 0.0) # if flanking_variant_count > 5, weight = 0.0

        if start_pos == 1: # read map to the left side
            if operate[-1] == 'M' or end_pos == ref_end_pos: # right side is open or till the end
                dict_allele_histogram[allele_name][start_pos-1:end_pos-1] += weight
        elif end_pos == ref_end_pos: # read map to the left side
            if operate[0] == 'M': # left side is open
                dict_allele_histogram[allele_name][start_pos-1:end_pos-1] += weight
        elif operate == ['M']: # read is including in the allele
            dict_allele_histogram[allele_name][start_pos-1:end_pos-1] += weight

        dict_allele_reads[allele_name].add(read_name)

def find_reference_allele(novel_allele_name, dict_ref_alleles):
    # known (Z49060|TRAV6*03|Homo), novel (|TRAV19*01/novel-0|)
    for ref_allele_name in dict_ref_alleles:
        ref_allele_pure_name = ref_allele_name.split('|')[1]
        # print(f'>>> novel_allele_name: {novel_allele_name}')
        if novel_allele_name == ref_allele_pure_name:
            # print(f'>>> ref_allele_name: {ref_allele_name}')
            return ref_allele_name
    # return None

def find_best_start_site(allele_sequence, reference_sequence, allele_name):
    best_start_site = 0
    position_score_dict={}

    for i in range(len(allele_sequence)-1):
        match_score = 0
        if i > (len(allele_sequence)-len(reference_sequence)):
            break
        sub_seq_of_allele_seq = allele_sequence[i:i+len(reference_sequence)]
        for j in range(len(sub_seq_of_allele_seq)-1):
            if sub_seq_of_allele_seq[j] == reference_sequence[j]:
                match_score += 1
        position_score_dict[i] = match_score

    if len(position_score_dict)==0: 
        print(f'Error: {allele_name} has no match score')
        print(f'>>> allele_sequence: {allele_sequence}')
        print(f'>>> reference_sequence: {reference_sequence}')    

    elif len(position_score_dict)>0:
        best_start_site = max(position_score_dict, key=position_score_dict.get)


    return best_start_site

def check_pair_end_coverage(allele_name, positions, primary_pair_dict, dict_allele_ori):
    """
    Check if pair-end reads cover the specified positions for a given allele.
    For single position: need at least one pair-end read covering it
    For two positions: need at least one pair-end read covering both
    For three+ positions: need at least one pair-end read covering adjacent position pairs (AB, BC, etc.)
    
    Args:
        allele_name: Name of the allele to check
        positions: List of positions to check
        primary_pair_dict: Dictionary of primary pair-end reads
        
    Returns:
        bool: True if the coverage criteria are met, False otherwise
    """
    if not positions:
        return True  # No positions to check
        
    # Convert positions to zero-based if they're not already (because python is zero-based)
    pos_zero_based = [p - 1 for p in positions]
    
    # Initialize a dictionary to track which positions are covered by which reads
    covered_positions = {pos: [] for pos in pos_zero_based}
    
    original_allele_name = find_reference_allele(allele_name, dict_allele_ori)
    for read_name, read_fields in primary_pair_dict.items():
        if len(read_fields) < 2:
            continue
            
        for field in read_fields:
            if field[2] != original_allele_name:
                continue
                
            start_pos = int(field[3]) - 1  # Convert to zero-based
            num_list, operate_list = parse_CIGAR(field[5])
            covered_length = sum(num for num, op in zip(num_list, operate_list) if op == 'M')
            end_pos = start_pos + covered_length
            
            # Check which positions are covered
            for pos in pos_zero_based:
                if start_pos <= pos < end_pos:
                    covered_positions[pos].append(read_name)
    
    # Print covered positions info
    for pos, reads in covered_positions.items():
        if reads:
            print(f"Position {pos + 1}: Covered by {', '.join(reads)}")
        else:
            print(f"Position {pos + 1}: 0")
    
    # If all positions were covered, return True
    return all(bool(reads) for reads in covered_positions.values())    
    

def add_read_to_perfect_list(list_perfect_fields, list_mismatch_fields, read_dict, FL_mark_dict):
    # list_perfect_fields: row in sam file
    # list_mismatch_fields: row in sam file
    # read_dict: key//read_name + tag; # value//ref_name, f'NM:i:{NM_tag}', f'MD:{MD_tag}', mismatches, ref_mismatch_pos, AS_tag
    # FL_mark_dict: key//ref_name; # value//FL_tag
    core_perfect_list = []
    salvaged_list = []
    discard_salvaged_list = []
    upper_limit_num_FLvariant = 5

    for row in list_mismatch_fields:
        field = row
        read_name = field[0]
        sum_flag = field[1]
        ref_name = field[2]
        read_name_unique = f'{read_name}_{sum_flag}_{ref_name}'
        tag_list = []
        
        if read_name_unique in read_dict:
            ref_name, NM_tag, MD_tag, mismatches, ref_mismatch_pos, AS_tag = read_dict[read_name_unique]
            core_range_start, core_range_end = FL_mark_dict[ref_name][0], FL_mark_dict[ref_name][1]

            core_variants = [1 if core_range_start <= pos <= core_range_end else 0 for pos in ref_mismatch_pos]
            flanking_variants = [1 if pos < core_range_start or pos > core_range_end else 0 for pos in ref_mismatch_pos]

            if sum(core_variants) == 0:
                core_perfect_list.append(f"{read_name_unique}_{ref_name}:{ref_mismatch_pos}->{core_variants}")

                # # Salvage the imperfect read with no core variants into the perfect list
                flanking_variant_count = sum(flanking_variants)
                if flanking_variant_count == 0:
                    list_perfect_fields.append(field)
                elif flanking_variant_count <= upper_limit_num_FLvariant:
                    salvaged_list.append((field, flanking_variant_count))
                else:
                    discard_salvaged_list.append((field, flanking_variant_count))

    # print(f'Num_of_CorePerfect_ImperectRead: {len(core_perfect_list)}')
    # print(f'Num_of_salvaged_read: {len(salvaged_list)}')
    # print(f'Num_of_discard_salvaged_read: {len(discard_salvaged_list)}')
    # print(f'>>> discard_salvaged_list: {discard_salvaged_list}')
    return salvaged_list 

def is_reclaimed_read(field, read_dict, FL_mark_dict):
    read_name = field[0]
    sum_flag = field[1]
    ref_name = field[2]
    read_name_unique = f'{read_name}_{sum_flag}_{ref_name}'

    if read_name_unique in read_dict:
        # read_dict[read_name] = [ref_name, f'NM:i:{NM_tag}', f'MD:{MD_tag}', mismatches, ref_mismatch_pos, AS_tag]
        _, _, _, _, ref_mismatch_pos, _ = read_dict[read_name_unique]
        if ref_name in FL_mark_dict:
            core_range_start, core_range_end = FL_mark_dict[ref_name]
            for variant_pos in ref_mismatch_pos:
                if core_range_start <= variant_pos <= core_range_end:
                    return False  # the read has a mismatch in the core region, so it is not reclaimed
            return True  # all mismatches are outside the core region, so it is reclaimed
        
    return False  # read not found in read_dict or no mismatches to check


if __name__ == '__main__':
    args = parse_args()
    fn_alleles = args.fn_alleles
    fn_sam = args.fn_sam
    thrsd = args.thrsd
    fo_calling_report  = args.fo_calling_report
    fo_grouping_pickle = args.fo_grouping_pickle
    fn_verify_annotation = args.fn_verify_annotation
    allele_name = args.allele_name
    fn_database_with_flanking = args.fn_database_with_flanking
    fn_database_without_flanking = args.fn_database_without_flanking
    fn_primary_aligned_pair = args.fn_primary_aligned_pair

    dict_allele, dict_novel_only_SEQ = parse_fasta3(fn_alleles) #include: novel_flanking (|TRAV19*01/novel-0|) + known_flanking (Z49060|TRAV6*03|Homo)
    dict_allele_ori = parse_fasta(fn_database_without_flanking) # known_w/o_flanking only
    FL_mark_dict = get_FL_mark_dict() # {'K02545|TRBD1*01|Homo': [200, 225], 'X02987|TRBD2*01|Homo': [200, 229]...}
    
    key_variants_to_check = positions_to_check()

    key_allele_dependent_variant = {'core':{}, 'flanking':{}, 'core_list':[], 'flanking_list':[]}

    
    differ_length_dict = {}

    # parse and add the novel allele into FL_mark_dict
    for novel_name, novel_seq in dict_novel_only_SEQ.items():
        novel_find_name = novel_name.split("|")[1][:-8]
        ref_name = find_reference_allele(novel_find_name, dict_allele_ori)
        if ref_name:
            ref_seq = dict_allele_ori.get(ref_name)
            best_start_site = find_best_start_site(novel_seq, ref_seq, novel_name)
            end_site = best_start_site + len(ref_seq) - 1
            FL_mark_dict[novel_name] = [best_start_site, end_site]

            # check the length of the novel allele and the reference allele
            reference_seq = dict_allele.get(ref_name)
            if len(novel_seq) != len(reference_seq):
                differ_length_dict[novel_name] = (f'{novel_name}:{len(novel_seq)}', f'{ref_name}:{len(reference_seq)}')
        else:
            eprint(f'Error: {novel_name} has no reference in the database_without_flanking')
            
    read_dict = parse_sam(fn_sam) # key//read_name; # value//ref_name, f'NM:i:{NM_tag}', f'MD:{MD_tag}', mismatches, ref_mismatch_pos
    dict_allele_histogram = { name.split()[0]:np.zeros(len(SEQ)) for name, SEQ in dict_allele.items() } # include: novel_flanking (|TRAV19*01/novel-0|) + known_flanking (Z49060|TRAV6*03|Homo)
    dict_allele_reads = { name.split()[0]:set() for name in dict_allele.keys() }
    print(f'>>> differ_length_dict: {differ_length_dict}')
  
    ################################
    ## parse the perfect sam file ##
    ################################
    list_perfect_fields, list_mismatch_fields = parse_perfect_sam_with_S(fn_sam)
    pri_list_perfect_fields, pri_list_mismatch_fields = parse_perfect_sam_with_S(fn_primary_aligned_pair)


    ### reclaim the reads that are in the mismatch list but should be in the perfect list
    salvaged_list = add_read_to_perfect_list(list_perfect_fields, list_mismatch_fields, read_dict, FL_mark_dict)
    histogram_read_depth(dict_allele_histogram, list_perfect_fields, dict_allele_reads, salvaged_list, thrsd)

    pri_salvaged_list = add_read_to_perfect_list(pri_list_perfect_fields, pri_list_mismatch_fields, read_dict, FL_mark_dict)
    histogram_read_depth(dict_allele_histogram, pri_list_perfect_fields, dict_allele_reads, pri_salvaged_list, thrsd, Pri=True)

    list_min_depth = []
    # list_for_check = []
    novel_list_for_check = []
    # print(f'>>> dict_allele_histogram_before:\n{dict_allele_histogram}')

    sample_list_for_check_FL_region = FL_region_for_check()
    for allele_name, histogram in dict_allele_histogram.items():
        if "novel" not in allele_name:
            allele_sequence = dict_allele.get(allele_name)
            reference_allele_name = allele_name
            reference_sequence = dict_allele_ori.get(reference_allele_name)
        elif "novel" in allele_name:
            novel_list_for_check.append(allele_name)
            novel_allele_name = allele_name.split("|")[1][:-8] # remove the "/novel-0" tag
            reference_allele_name = find_reference_allele(novel_allele_name, dict_allele_ori)
            if reference_allele_name:
                reference_sequence = dict_allele_ori.get(reference_allele_name)
            else: print(f'Error: {allele_name} has no reference in the database_without_flanking')
        else: print("ERROR while checking allele_name")

        best_start_site = find_best_start_site(allele_sequence, reference_sequence, allele_name)
        end_site = best_start_site + len(reference_sequence) #- 1
        if end_site > len(histogram):
            min_depth = 0
            print(f'ERROR: {allele_name} has end_site > len(histogram)')
        elif allele_name in sample_list_for_check_FL_region:
            min_depth = min(histogram[best_start_site:len(histogram)])

        else:
            min_depth = min(histogram[best_start_site:end_site])


        list_min_depth.append([min_depth, None, allele_name.split('|')[1]])


    ### check key variants
    # Parse primary pair-end reads
    primary_pair_dict = parse_pair_sam(fn_primary_aligned_pair)
    # Check pair-end coverage for specified alleles
    updated_list_min_depth = []
    for element in list_min_depth:
        min_depth, tag, name = element
        
        # Check if this allele needs the pair-end check
        if min_depth != 0:
            for check_name, positions in key_variants_to_check.items():
                if check_name == name:  # Check if the allele name matches or contains the check name
                    # Perform the pair-end coverage check
                    if not check_pair_end_coverage(name, positions, primary_pair_dict, dict_allele_ori):
                        print(f">>> Allele {name} failed pair-end coverage check. Original depth: {min_depth}, Setting to 0.1111")
                        min_depth = 0.1111

                    break
        
        updated_list_min_depth.append([min_depth, tag, name])
    
    # Replace the original list with the updated one
    list_min_depth = updated_list_min_depth


    set_annotation = set()
    if fn_verify_annotation:
        with open(fn_verify_annotation, "r") as f_v:
            for line in f_v:
                set_annotation.add(line.split(',')[0])

    for idx, ele in enumerate(list_min_depth):
        min_depth, tag, name = ele
        if "corrected" in name:
            list_min_depth[idx][1] = 4
            continue
        if name in set_annotation:
            if min_depth > 0.5:
                list_min_depth[idx][1] = 3
            else:
                list_min_depth[idx][1] = 1
        else:
            if min_depth > 0.5:
                list_min_depth[idx][1] = 2
            else:
                list_min_depth[idx][1] = 0

    f_oc = open(fo_calling_report,'w')
    for element in sorted(list_min_depth, reverse=True):
        min_depth, tag, name = element
        if fn_verify_annotation:
            if tag==0:
                continue
            f_oc.write(name + '\t' + str(min_depth) + '\t' + str(tag))
        else:
            f_oc.write(name + '\t' + str(min_depth))
        f_oc.write('\n')
    f_oc.close()

    f_op = open(fo_grouping_pickle, 'wb')
    pickle.dump(dict_allele_reads, f_op)

    

