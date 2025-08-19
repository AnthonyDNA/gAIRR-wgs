import argparse
import pickle
import os
import numpy as np
import sys
from test2 import parse_sam, get_FL_mark_dict
from test3 import parse_fasta2, parse_fasta3
from collections import defaultdict
# make sure the package modules is in the path
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.dirname(__file__)+'/..')

from utils import get_reverse_complement
from filter_corrected_alleles import parse_perfect_sam, parse_fasta
from parse_contig_realign import parse_CIGAR

import re
def extract_mismatches(md_tag):
    # 使用正規表達式找出數字和字母的組合
    matches = re.findall(r'(\d+)([A-Z]+)', md_tag)
    
    mismatch_positions = []
    current_position = 0
    
    # 遍歷找到的組合，計算 mismatch 的位置並添加到列表中
    for num, letters in matches:
        current_position += int(num)
        mismatch_positions.append(current_position + 1)  # 加1表示 mismatch 的位置
        current_position += len(letters)
    
    return mismatch_positions
def parse_sam(sam_file):
    read_dict = {} #key: read name; value: read seq, ref name, ref relative pos (ref pos+MD_pos), CIGAR_mid, MD_pos
    with open(sam_file, 'r') as f:
        for line in f:
            if not line.startswith('@'):
                sumflag = line.strip().split()[1]
                read_seq = line.strip().split()[9]
                ref_name = line.strip().split()[2]
                read_name = f'{line.strip().split()[0]}_{sumflag}_{ref_name}'
                ref_pos = line.strip().split()[3]
                cigar = line.strip().split()[5]
                mismatches=0
                ref_mismatch_pos=0
                AS_tag = 0
                NM_tag = None
                MD_tag = None
                if line.strip().split()[11] != "" or line.strip().split()[12] != "":
                    if "NM" in line.strip().split()[11]:
                        NM_tag = line.strip().split()[11].split(':')[-1]
                        if "MD" in line.strip().split()[12]:
                            MD_tag = line.strip().split()[12].split(':')[-1]
                            mismatches = extract_mismatches(MD_tag)
                            ref_mismatch_pos = [int(ref_pos) + int(mismatch) -2 for mismatch in mismatches]  #python seen position 1 as 0
                            # read_dict[read_name] = [ref_name, f'NM:i:{NM_tag}', f'MD:{MD_tag}', mismatches, ref_mismatch_pos]
                    if "AS" in line.strip().split():
                        AS_idx = line.strip().split().index("AS:i:")
                        AS_tag = line.strip().split()[AS_idx].split(':')[-1]
                        AS_tag = int(AS_tag)
                read_dict[read_name] = [ref_name, f'NM:i:{NM_tag}', f'MD:{MD_tag}', mismatches, ref_mismatch_pos, AS_tag]

    return(read_dict)

def parse_fasta3(fn_fasta):
    '''parse the fasta file into a dictionary'''
    # dict_name_SEQ {}
    #  - keys: seq_name
    #  - values: seq_SEQ
    dict_name_SEQ = {}
    dict_novel_only_SEQ = {}
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
        
        # add novel into dict_novel_only_SEQ
        for seq_name in dict_name_SEQ:
            if "novel" in seq_name:
                dict_novel_only_SEQ[seq_name] = dict_name_SEQ[seq_name]

    return dict_name_SEQ, dict_novel_only_SEQ

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
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_perfect_sam_with_S(fn_sam, read_dict):
    list_perfect_fields = []
    list_mismatch_fields = []
    pair_end_dict = {}
    # MAPQ_dict = {}
    with open(fn_sam, 'r') as f_s:
        for line in f_s:
            if line[0] != '@':
                fields = line.strip().split()
                # if "NM:i:0" in line or "NM:i:1" in line or "NM:i:2" in line:
                # if "NM:i:0" in line:
                if "NM:i:0" in line:
                    list_perfect_fields.append(fields)
                else:
                    list_mismatch_fields.append(fields)
                read_name = fields[0]
                sum_flag = fields[1]
                ref_name = fields[2]
                read_unique_name = f'{read_name}_{sum_flag}_{ref_name}'
                AS_tag = read_dict[read_unique_name][5]
                # mapQ = fields[4]
                # if mapQ not in MAPQ_dict:
                #     MAPQ_dict[mapQ] = 0
                # MAPQ_dict[mapQ] += 1
                # cigar = fields[5]
                
                if AS_tag == None:
                    AS_tag = 0
                fields.append(AS_tag)
                if read_name not in pair_end_dict:
                    pair_end_dict[read_name]=[]

                pair_end_dict[read_name].append(tuple(fields)) # fields is a list, so need to convert to tuple
                # in pair_end_dict, the last element is AS_tag
                
                

    # print(f'MAPQ: {MAPQ_dict}')
    print(f'Total Num of Cluster: {len(pair_end_dict)}')            

    # Parse the pair_end_dict
    len_count_dict = defaultdict(int)
    count_dict = defaultdict(lambda: defaultdict(int))
    for read_name, values in pair_end_dict.items(): # each read has multiple alignments
        value_len = len(values) # the number of alignment sites
        len_count_dict[value_len] += 1

        ref_names = [t[2] for t in values]

        ref_set_len = len(set(ref_names)) # Num of types of ref for each read
        ref_count = defaultdict(int)
        for ref in ref_names:
            ref_count[ref] += 1
        
        ref_count_sorted = tuple(sorted(ref_count.values(), reverse=True))
        count_dict[value_len][ref_count_sorted] += 1
    
    print(f'Num of each cluster: {len_count_dict}')
    # print(count_dict)



    return list_perfect_fields, list_mismatch_fields, pair_end_dict

def process_pair(dict_allele_histogram, pair1, pair2):
    allele_name1, allele_name2 = pair1[2], pair2[2]
    if allele_name1 != allele_name2:
        return True
    else:
        # parse_CIGAR: num_list, operate_list
        num_list1, operate_list1 = parse_CIGAR(pair1[5])
        num_list2, operate_list2 = parse_CIGAR(pair2[5])
        start_pos1 = int(pair1[3])
        end_pos1 = start_pos1 + sum(num for num, op in zip(num_list1, operate_list1) if op == 'M')
        start_pos2 = int(pair2[3])
        end_pos2 = start_pos2 + sum(num for num, op in zip(num_list2, operate_list2) if op == 'M')

        # start_pos1, end_pos1 = int(pair1[3]), int(pair1[3]) + sum([x[0] for x in parse_CIGAR(pair1[5])[0] if x[1] == 'M']) # ideally, the end_pos should minus 1
        # start_pos2, end_pos2 = int(pair2[3]), int(pair2[3]) + sum([x[0] for x in parse_CIGAR(pair2[5])[0] if x[1] == 'M'])
        coverage_range = range(min(start_pos1, start_pos2), max(end_pos1, end_pos2)) # however, the end_pos should plus 1 to make sure the end pos is included
        dict_allele_histogram[allele_name1][list(coverage_range)] += 1
        return False
    # except:
    #     # eprint("SKIP", pair1, pair2)
    #     return True

def histogram_read_depth(dict_allele_histogram, list_perfect_fields, dict_allele_reads, pair_end_dict, read_dict, thrsd=100, min_AS=30):
    for fields in list_perfect_fields: # fields: fields = line.strip().split()
        allele_name = fields[2]
        read_name = fields[0]
        flag = fields[1]
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
    bonus_read_list = []
    for read_name, read_fields in pair_end_dict.items():
        if len(read_fields) == 1:
            continue # skip the single-end reads

        # read_fields in pair_end_dict is a list, each element (read_field) in read_fields is a tuple, one element means that one read align to one position
        # AS_tag is the last element in the tuple, so the index is -1; ref_name is the 3rd element, so the index is 2
        
        pairs_with_AS = [(read_field, int(read_field[-1]), read_field[2], read_field[8]) for read_field in read_fields]
        pairs_with_AS.sort(key=lambda x: x[1], reverse=True)
        if len(read_fields) == 2: # read align to two positions, may be forward and reverse
            field1, field2 = [p[0] for p in pairs_with_AS] # to get the read field (tuple), which is the first element in the tuple
            if field1[2] == field2[2]: # if the two reads align to the same reference
                tlen_1, tlen_2 = int(pairs_with_AS[0][3]), int(pairs_with_AS[1][3])
                as_1, as_2 = pairs_with_AS[0][1], pairs_with_AS[1][1]
                # using tlen to determine the forward and reverse
                if (tlen_1 * tlen_2) < 0 and as_1 >= min_AS and as_2 >= min_AS:
                    process = process_pair(dict_allele_histogram, field1, field2)
                    if not process:
                        data_array = (field1[2], f'{field1[0]}_{field1[1]}', f'{field2[0]}_{field2[1]}')
                        bonus_read_list.append(data_array)
            # # if AS > min_AS, then process the pair
            # if pairs_with_AS[0][1] >= min_AS and pairs_with_AS[1][1] >= min_AS: 
            #     process = process_pair(dict_allele_histogram, field1, field2)
        elif len(read_fields) > 2: # read align to more than two positions
            high_AS_fields = [p for p in pairs_with_AS if p[1] >= min_AS] # p (tuple) is each element in pairs_with_AS, p[1] is AS_tag; high_AS_fields is a list of tuples

            while high_AS_fields: # parse each alignment position (tuple); tuple: (read_field, AS_tag, ref_name, TLEN_num)
                best_pair = max(high_AS_fields, key=lambda x: (x[1], pairs_with_AS.index(x))) # get the best pair with the highest AS_tag; best_pair is a tuple (read info) 
                ref_name = best_pair[0][2]
                tlen_num = int(best_pair[3])
                mate_pairs = [p for p in high_AS_fields if p[0][2] == ref_name and p != best_pair and (tlen_num * int(p[3])) < 0] # "mate_pairs" is a list of tuples, each tuple is paired read aligning to the same reference
                if mate_pairs:
                    process = True
                    for mate_pair in mate_pairs:
                        if not process: # process is False, meaning the mate_pair was successfully found
                            break
                        else:
                            process == process_pair(dict_allele_histogram, best_pair[0], mate_pair[0]) #best_pair[0] and mate_pair[0] are "read field (tuple format)""
                            if not process: # process is False, meaning the mate_pair was successfully found; then remove the best_pair and mate_pair from high_AS_fields
                                high_AS_fields.remove(best_pair)
                                high_AS_fields.remove(mate_pair)
                                data_array = (best_pair[2], f'{best_pair[0][0]}_{best_pair[0][1]}', f'{mate_pair[0][0]}_{mate_pair[0][1]}')
                                bonus_read_list.append(data_array)
                                break
                else:
                    high_AS_fields.remove(best_pair)
    print(f'>>> bonus_read_list: {bonus_read_list}')
            # else:
            #     best_pair = pairs_with_mapq[0]
            #     ref_name = best_pair[0][2]
            #     mate_pair = [p for p in pairs_with_mapq[1:] if p[0][2] == ref_name]
            #     for field in mate_pair:
            #         process_pair(dict_allele_histogram, best_pair[0], field[0])

def find_reference_allele(novel_allele_name, dict_ref_alleles):
    # known (Z49060|TRAV6*03|Homo), novel (|TRAV19*01/novel-0|)
    for ref_allele_name in dict_ref_alleles:
        ref_allele_pure_name = ref_allele_name.split('|')[1]
        # print(f'>>> novel_allele_name: {novel_allele_name}')
        if novel_allele_name == ref_allele_pure_name:
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
    if len(position_score_dict)==0: print(f'Error: {allele_name} has no match score')
    elif len(position_score_dict)>0:
        best_start_site = max(position_score_dict, key=position_score_dict.get)
    return best_start_site
def add_read_to_perfect_list(list_perfect_fields, list_mismatch_fields, read_dict, FL_mark_dict):
    # list_perfect_fields: row in sam file
    # list_mismatch_fields: row in sam file
    # read_dict: key//read_name + tag; # value//ref_name, f'NM:i:{NM_tag}', f'MD:{MD_tag}', mismatches, ref_mismatch_pos, AS_tag
    # FL_mark_dict: key//ref_name; # value//FL_tag
    reclaim_list = []
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
            for variant_pos in ref_mismatch_pos:
                if core_range_start <= variant_pos <= core_range_end:
                    tag_list.append(1)
                else:
                    tag_list.append(0)
            if sum(tag_list) == 0:
                reclaim_list.append(f"{read_name_unique}_{ref_name}:{ref_mismatch_pos}->{tag_list}")
                list_perfect_fields.append(field)
    # print(f'>>> Reclaim_list: {reclaim_list}')
    print(f'Num_of_reclaimed_read: {len(reclaim_list)}')

if __name__ == '__main__':
    args = parse_args()
    fn_alleles = args.fn_alleles #material_allele_fasta
    fn_sam = args.fn_sam #sam file
    thrsd = args.thrsd
    fo_calling_report  = args.fo_calling_report
    fo_grouping_pickle = args.fo_grouping_pickle
    fn_verify_annotation = args.fn_verify_annotation
    allele_name = args.allele_name
    fn_database_with_flanking = args.fn_database_with_flanking #database_with_flanking.fasta (parsed_material)
    fn_database_without_flanking = args.fn_database_without_flanking #database_without_flanking.fasta (raw_material)

    dict_allele, dict_novel_only_SEQ = parse_fasta3(fn_alleles) #include: novel_flanking (|TRAV19*01/novel-0|) + known_flanking (Z49060|TRAV6*03|Homo)

    dict_allele_ori = parse_fasta(fn_database_without_flanking) # known_w/o_flanking only
    FL_mark_dict = get_FL_mark_dict() # {'K02545|TRBD1*01|Homo': [200, 225], 'X02987|TRBD2*01|Homo': [200, 229]...}
    # parse and add the novel allele into FL_mark_dict
    for novel_name, novel_seq in dict_novel_only_SEQ.items():
        novel_find_name = novel_name.split("|")[1][:-8]
        ref_name = find_reference_allele(novel_find_name, dict_allele_ori)
        if ref_name:
            ref_seq = dict_allele_ori.get(ref_name)
            best_start_site = find_best_start_site(novel_seq, ref_seq, novel_name)
            end_site = best_start_site + len(ref_seq) - 1
            FL_mark_dict[novel_name] = [best_start_site, end_site]
        else:
            eprint(f'Error: {novel_name} has no reference in the database_without_flanking')
    read_dict = parse_sam(fn_sam) # key//read_name; # value//ref_name, f'NM:i:{NM_tag}', f'MD:{MD_tag}', mismatches, ref_mismatch_pos

    # print(f'>>> dict_allele_ori:\n{dict_allele_ori}')
    # '''collect the reference sequence for novel alleles'''
    # ## collect the novel allele name from
    # novel_allele_names = [name for name in dict_allele if "novel" in name]
    # ## collect the reference sequence for novel alleles through the novel allele name
    # dict_novel_project_to_ref = {}
    # for novel_allele_name in novel_allele_names:
    #     reference_allele_name = find_reference_allele(novel_allele_name, dict_allele_ori)
    #     if reference_allele_name:
    #         reference_sequence = dict_allele_ori.get(reference_allele_name)
    #         dict_novel_ref_sequences[novel_allele_name] = reference_sequence
    #     else:
    #         eprint(f'Error: {novel_allele_name} has no reference in the database_without_flanking')
    # dict_allele_merge = {**dict_allele_ori, **dict_novel_project_to_ref}

    dict_allele_histogram = { name.split()[0]:np.zeros(len(SEQ)) for name, SEQ in dict_allele.items() }
    # include: novel_flanking (|TRAV19*01/novel-0|) + known_flanking (Z49060|TRAV6*03|Homo)

    dict_allele_reads = { name.split()[0]:set() for name in dict_allele.keys() }

    list_perfect_fields, list_mismatch_fields, pair_end_dict = parse_perfect_sam_with_S(fn_sam, read_dict)
    add_read_to_perfect_list(list_perfect_fields, list_mismatch_fields, read_dict, FL_mark_dict)
    histogram_read_depth(dict_allele_histogram, list_perfect_fields, dict_allele_reads, pair_end_dict, read_dict, thrsd, min_AS=30)
    # print(f'>>> dict_allele_reads:\n{dict_allele_reads}')
    # print(f'>>> dict_allele_histogram:\n{dict_allele_histogram}')

    # for allele_name, histogram in dict_allele_histogram.items():
    #     allele_sequence = dict_allele.get(allele_name)
    #     reference_allele_name = find_reference_allele(allele_name, dict_allele_ori)
    #     if reference_allele_name:
    #         reference_sequence = dict_allele_ori.get(reference_allele_name)
    #     else: print(f'Error: {allele_name} has no reference in the database_without_flanking')

    # list_min_depth = [ [min(histogram), None, name.split('|')[1]] for name, histogram in  dict_allele_histogram.items() ]
    list_min_depth = []
    list_for_check = []
    novel_list_for_check = []
    # print(f'>>> dict_allele_histogram_before:\n{dict_allele_histogram}')
    for allele_name, histogram in dict_allele_histogram.items():
        # print("checking allele_name: ", allele_name) ## allele_name: known (Z49060|TRAV6*03|Homo), novel (|TRAV19*01/novel-0|)
        list_for_check.append(allele_name)
        if "novel" not in allele_name:
            allele_sequence = dict_allele.get(allele_name)
            reference_allele_name = allele_name
            reference_sequence = dict_allele_ori.get(reference_allele_name)
        elif "novel" in allele_name:
            novel_list_for_check.append(allele_name)
            novel_allele_name = allele_name.split("|")[1][:-8]
            reference_allele_name = find_reference_allele(novel_allele_name, dict_allele_ori)
            if reference_allele_name:
                reference_sequence = dict_allele_ori.get(reference_allele_name)
                # print(f'Parsed Ref of novel allele: {allele_name} -> {reference_allele_name}')
            else: print(f'Error: {allele_name} has no reference in the database_without_flanking')
        else: print("ERROR while checking allele_name")
        best_start_site = find_best_start_site(allele_sequence, reference_sequence, allele_name)
        end_site = best_start_site + len(reference_sequence) - 1
        min_depth = min(histogram[best_start_site:end_site+1])
        list_min_depth.append([min_depth, None, allele_name.split('|')[1]])
    # print(f'>>> list_for_check: {list_for_check}')
    # print(f'>>> novel_list_for_check: {novel_list_for_check}')

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
            if min_depth > 0:
                list_min_depth[idx][1] = 3
            else:
                list_min_depth[idx][1] = 1
        else:
            if min_depth > 0:
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

    

