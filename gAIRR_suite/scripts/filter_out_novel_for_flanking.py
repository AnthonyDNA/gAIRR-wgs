import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fcafil', '--fn_corrected_alleles_raw', 
        help='input corrected_alleles_raw.fasta'
    )
    parser.add_argument(
        '-fdbfnk', '--fn_database_with_flanking', 
        help='input database_with_flanking.fasta'
    )
    parser.add_argument(
        '-fdbraw', '--fn_database_without_flanking', 
        help='input raw_database.fasta'
    )
    return parser.parse_args()

def parse_fasta2(fn_fasta):
    ''' parse the fasta file into a dictionary
    # dict_name_SEQ {}
    #  - keys: seq_name -> X72718|TRBV20/OR9-2*02|Homo
    #  - values: seq_SEQ -> all small letters
    '''
    dict_name_SEQ = {}
    with open(fn_fasta, 'r') as f_f:
        seq_name = ""
        seq_SEQ = ""
        for line in f_f:
            if line[0] == '>':
                seq_name = line.strip()[1:]
                try:
                    seq_name = seq_name.split()[0]
                except:
                    pass
                seq_SEQ = ""
            else:
                seq_SEQ += line.strip()
            dict_name_SEQ[seq_name] = seq_SEQ

    return dict_name_SEQ

def parse_fasta3(fn_fasta):
    ''' parse the fasta file into a dictionary
    # novel_allele_dict {}
    #  - keys: novel_seq_name -> X72718|TRBV20/OR9-2*02|Homo/novel1
    #  - values: novel_seq -> contain Capital letter, marking the variant
    '''
    novel_allele_dict = {}
    novel_allele_name_list = []
    with open(fn_fasta, 'r') as f_f:
        seq_name = ""
        novel_seq = ""
        for line in f_f:
            if line[0] == '>':
                seq_name = line.strip()[1:]
                try:
                    seq_name = seq_name.split()[0]
                except:
                    pass
                novel_seq = ""
            else:
                novel_seq += line.strip()

            if "novel" in seq_name:
                novel_allele_dict[seq_name] = novel_seq
                novel_allele_name_list.append(seq_name)
    return novel_allele_dict, novel_allele_name_list

def extract_middle_part(seq_name):
    """Extract the middle part of sequence name after splitting by '|'"""
    parts = seq_name.split('|')
    if len(parts) >= 2:
        return parts[1]
    return seq_name

def find_matching_database_name(core_seq, db_without_flanking_dict):
    """Find the database name that matches the core sequence"""
    for db_name, db_seq in db_without_flanking_dict.items():
        if core_seq.lower() == db_seq.lower():
            return db_name
    return None

def check_if_name_exists_in_corrected_file(target_middle_part, novel_allele_name_list):
    """Check if the target middle part exists in the corrected alleles file"""
    for novel_name in novel_allele_name_list:
        novel_middle_part = extract_middle_part(novel_name)
        if target_middle_part == novel_middle_part:
            return True, novel_name
    return False, None

def record_differences(sequence1, sequence2):
    sequence1=sequence1.lower()
    sequence2=sequence2.lower()
    differences = []
    for i, (char1, char2) in enumerate(zip(sequence1, sequence2)):
        if char1 != char2:
            differences.append(i)
    return differences

def record_differences_aligned(sequence1, sequence2, best_start_site, best_start_site_on_other):
    differences = []
    if len(sequence1) > len(sequence2):
        sub_seq1 = sequence1[best_start_site:(best_start_site+len(sequence2))]
        sub_seq2 = sequence2
        if len(sub_seq1) != len(sub_seq2): 
            print("RecordError//seq1>seq2: the length of sub_seq1 is not the same as sub_seq2")
        for j, (char1, char2) in enumerate(zip(sub_seq1, sub_seq2)):
            if char1 != char2:
                differences.append(best_start_site+j)
    elif len(sequence2) > len(sequence1):
        sub_seq1 = sequence1
        sub_seq2 = sequence2[best_start_site_on_other:(best_start_site_on_other+len(sequence1))]
        if len(sub_seq1) != len(sub_seq2): 
            print("RecordError//seq2>seq1: the length of sub_seq1 is not the same as sub_seq2")
        for j, (char1, char2) in enumerate(zip(sub_seq1, sub_seq2)):
            if char1 != char2:
                differences.append(j)
    elif len(sequence1) == len(sequence2):
        sub_seq1 = sequence1[best_start_site:]
        sub_seq2 = sequence2[:(len(sequence2)-best_start_site)]
        if len(sub_seq1) != len(sub_seq2): 
            print("RecordError//seq1=seq2: the length of sub_seq1 is not the same as sub_seq2")
        for j, (char1, char2) in enumerate(zip(sub_seq1, sub_seq2)):
            if char1 != char2:
                differences.append(best_start_site+j)
    return differences

def find_best_start_site(sequence1, sequence2):
    ### find the best match between the two sequences
    best_start_site = 0
    best_start_site_on_other = 0
    position_score_dict={}
    position_score_on_other_dict={}
    if len(sequence1) > len(sequence2):
        for i in range(len(sequence1)-1):
            match_score = 0
            if i > (len(sequence1)-len(sequence2)):
                break
            sub_seq1=sequence1[i:(i+len(sequence2))]
            if len(sub_seq1) != len(sequence2):
                print("sequence1 > sequence2//Error: the length of sub_seq1 is not the same as sequence2")
            for j in range(len(sub_seq1)-1):
                if sub_seq1[j] == sequence2[j]:
                    match_score += 1
            position_score_dict[i] = match_score
    elif len(sequence2) > len(sequence1):
        for i in range(len(sequence2)-1):
            match_score = 0
            if i > (len(sequence2)-len(sequence1)):
                break
            sub_seq2=sequence2[i:(i+len(sequence1))]
            if len(sub_seq2) != len(sequence1):
                print("sequence2 > sequence1//Error: the length of sub_seq2 is not the same as sequence1")
            for j in range(len(sub_seq2)-1):
                if sub_seq2[j] == sequence1[j]:
                    match_score += 1
            position_score_on_other_dict[i] = match_score
    elif len(sequence1) == len(sequence2):
        for i in range(len(sequence1)-1):
            match_score = 0
            sub_seq1 = sequence1[i:]
            sub_seq2 = sequence2[:(len(sequence2)-i)]
            if len(sub_seq1) != len(sub_seq2): 
                print("sequence1 = sequence2//Error: the length of sub_seq1 is not the same as sub_seq2")
            for j in range(len(sub_seq1)-1):
                if sub_seq1[j] == sub_seq2[j]:
                    match_score += 1
            position_score_dict[i] = match_score
    if len(position_score_dict) == 0 and len(position_score_on_other_dict) == 0:
        print("ERROR: position_score_dict or on_other is empty")
    elif len(position_score_dict) > 0:
        best_start_site = max(position_score_dict, key=position_score_dict.get)
    elif len(position_score_on_other_dict) > 0:
        best_start_site_on_other = max(position_score_on_other_dict, key=position_score_on_other_dict.get)
        
    return best_start_site, best_start_site_on_other

def find_variant_btw_seq_with_flanking(novel_allele_seq, database_ref_allele_seq_with_flanking, novel_allele_name, database_ref_allele_name):
    '''find the variant position between novel and database_with_flanking (ideally they have the same seq length)'''
    if len(novel_allele_seq) != len(database_ref_allele_seq_with_flanking):
        print(f'CAUTION: {novel_allele_name}:{len(novel_allele_seq)} and {database_ref_allele_name}:{len(database_ref_allele_seq_with_flanking)} have different length')
        best_start_site, best_start_site_on_other = find_best_start_site(novel_allele_seq, database_ref_allele_seq_with_flanking)
        variant_pos_list = record_differences_aligned(novel_allele_seq, database_ref_allele_seq_with_flanking, best_start_site, best_start_site_on_other)
    else:
        variant_pos_list = record_differences(novel_allele_seq, database_ref_allele_seq_with_flanking)
    return variant_pos_list
      
def check_variant_pos(novel_allele_dict, novel_allele_name_list, db_without_flanking_dict, db_with_flanking_dict):
    novel_alleles_for_retention = []
    discard_list = []
    core_seq_dict = {}
    core_seq_of_collected_novel_allele = {}
    
    '''collect the core sequence of the novel allele for later comparison'''
    for novel_allele_name in novel_allele_name_list:
        novel_allele_seq = novel_allele_dict[novel_allele_name]
        database_ref_allele_name = novel_allele_name[:-7]
        database_ref_allele_seq_with_flanking = db_with_flanking_dict[database_ref_allele_name]
        index = database_ref_allele_seq_with_flanking.find(db_without_flanking_dict[database_ref_allele_name])
        core_seq = novel_allele_seq[index:index+len(db_without_flanking_dict[database_ref_allele_name])]
        core_seq_dict[novel_allele_name] = core_seq
    
    '''check novel alleles iteratively by comparing with the database_with_flanking and database_without_flanking'''
    for novel_allele_name in novel_allele_name_list:
        novel_allele_seq = novel_allele_dict[novel_allele_name]
        database_ref_allele_name = novel_allele_name[:-7]
        database_ref_allele_seq_with_flanking = db_with_flanking_dict[database_ref_allele_name]
        database_ref_allele_seq_without_flanking = db_without_flanking_dict[database_ref_allele_name]
        
        '''check each variant whether it is in the alignment position or not'''
        variant_pos_list = find_variant_btw_seq_with_flanking(novel_allele_seq, database_ref_allele_seq_with_flanking, novel_allele_name, database_ref_allele_name)
        
        index = database_ref_allele_seq_with_flanking.find(database_ref_allele_seq_without_flanking)
        start_pos = index
        end_pos = index+len(database_ref_allele_seq_without_flanking)-1
        
        core_seq = core_seq_dict[novel_allele_name]
        
        annotate_list = []
        annotate_list_for_print = []
        
        '''deal with the same core sequence with its ref (i.e. False Positive)'''
        if core_seq.lower() == database_ref_allele_seq_without_flanking.lower():
            discard_list.append(novel_allele_name)
            print(f'\t\tDISCARD {novel_allele_name} due to the same core sequence with REF:{database_ref_allele_name}')
            print(f'\t\t\tQ: {core_seq}\n\t\t\tR: {database_ref_allele_seq_without_flanking}')
        elif any(core_seq.lower() == db_seq.lower() for db_seq in db_without_flanking_dict.values()):
            # find matching database name
            matching_db_name = find_matching_database_name(core_seq, db_without_flanking_dict)
            if matching_db_name:
                # extract the middle part of the matching database name for comparison
                matching_middle_part = extract_middle_part(matching_db_name)
                exists_in_corrected, existing_novel_name = check_if_name_exists_in_corrected_file(matching_middle_part, novel_allele_name_list)
                
                if exists_in_corrected:
                    # if the corresponding name exists in the corrected file, discard the novel allele
                    discard_list.append(novel_allele_name)
                    print(f'\t\tDISCARD {novel_allele_name} due to the same core sequence with database entry: {matching_db_name}')
                    print(f'\t\t\tQ: {core_seq}\n\t\t\tR: {db_without_flanking_dict[matching_db_name]}')
                    print(f'\t\t\tNOTE: Corresponding name "{matching_middle_part}" already EXISTS in corrected file as: {existing_novel_name}')
                else:
                    # if the corresponding name does not exist in the corrected file, retain the novel allele with the complete database name
                    print(f'\t\tRENAME & RETAIN {novel_allele_name} -> {matching_db_name} (same core sequence with database entry)')
                    print(f'\t\t\tQ: {core_seq}\n\t\t\tR: {db_without_flanking_dict[matching_db_name]}')
                    print(f'\t\t\tNOTE: Corresponding name "{matching_middle_part}" does NOT exist in corrected file, so retaining with complete database name: {matching_db_name}')
                    
                    # prepare the output line with the new name
                    name_for_output = f'{matching_db_name}_DM_'
                    output_in_line = f'>{name_for_output}\n{novel_allele_dict[novel_allele_name]}\n'
                    # check if the output line is already in the retention list
                    existing_names = [line.split('\n')[0][1:] for line in novel_alleles_for_retention if line.startswith('>')]
                    if name_for_output not in existing_names:  # 修正：檢查完整的輸出名稱
                        novel_alleles_for_retention.append(output_in_line)
                        core_seq_of_collected_novel_allele[name_for_output] = core_seq  # 修正：使用輸出名稱作為key

            else:
                discard_list.append(novel_allele_name)
                print(f'\t\tDISCARD {novel_allele_name} due to the same core sequence with other database entries (cannot find specific match)')
                print(f'\t\t\tQ: {core_seq}')
        else:
            '''annotate the variant position'''
            for variant_pos in variant_pos_list:
                if variant_pos < start_pos or variant_pos > end_pos:
                    line = f'{variant_pos}:0'
                    annotate_list_for_print.append(line)
                    annotate_list.append("0")
                else:
                    line = f'{variant_pos}:1'
                    annotate_list_for_print.append(line)
                    annotate_list.append("1")
            
            if "1" in annotate_list and novel_allele_name not in discard_list:
                if novel_allele_name in discard_list:
                    print(f'\t\tDISCARD {novel_allele_name} due to already discarded')
                elif novel_allele_name in novel_alleles_for_retention:
                    print(f'\t\tDISCARD {novel_allele_name} due to already collected')
                elif novel_allele_name in discard_list and novel_allele_name in novel_alleles_for_retention:
                    print(f'\t\tcontradictoryERROR: {novel_allele_name} is in both discard_list and novel_alleles_for_retention')
                elif novel_allele_name not in discard_list:
                    if len(novel_alleles_for_retention) == 0:
                        output_in_line = f'>{novel_allele_name}\n{novel_allele_dict[novel_allele_name]}\n'
                        if output_in_line not in novel_alleles_for_retention and novel_allele_name not in discard_list:
                            novel_alleles_for_retention.append(output_in_line)
                            core_seq_of_collected_novel_allele[novel_allele_name] = core_seq
                            print(f'\t\tPrimaryRETAIN {novel_allele_name}: {annotate_list_for_print}')
                    else:
                        retain_mark = True
                        core_seq_of_collected_novel_allele_copy = core_seq_of_collected_novel_allele.copy()
                        for allele_N, core_S in core_seq_of_collected_novel_allele_copy.items():
                            if novel_allele_name == allele_N:
                                print(f'CAUTION: {novel_allele_name} has been already collected in the novel_alleles_for_retention list')
                                print(f'\t\tChecking the core sequence...')
                                if core_S.lower() != core_seq.lower():
                                    print(f'\t\tcore_seqERROR: {novel_allele_name} has the same name but different core sequence; corrected_raw.fasta should be checked')
                            elif novel_allele_name != allele_N:
                                if core_S.lower() == core_seq.lower():
                                    discard_list.append(novel_allele_name)
                                    print(f'\t\tDISCARD {novel_allele_name} due to the same core sequence with {allele_N}, which has been already collected')
                                    print(f'\t\t\tQ: {core_seq}\n\t\t\tR: {core_S}')
                                    retain_mark = False
                                    break
                        if retain_mark:
                            output_in_line = f'>{novel_allele_name}\n{novel_allele_dict[novel_allele_name]}\n'
                            if output_in_line not in novel_alleles_for_retention and novel_allele_name not in discard_list:
                                novel_alleles_for_retention.append(output_in_line)
                                core_seq_of_collected_novel_allele[novel_allele_name] = core_seq
                                print(f'\t\tSecondaryRETAIN {novel_allele_name}: {annotate_list_for_print}')
        
    return novel_alleles_for_retention

if __name__ == "__main__":
    args = parse_args()
    corrected_alleles_raw_fasta = args.fn_corrected_alleles_raw
    database_with_flanking_fasta = args.fn_database_with_flanking
    database_without_flanking_fasta = args.fn_database_without_flanking

    db_without_flanking_dict = parse_fasta2(database_without_flanking_fasta)
    db_with_flanking_dict = parse_fasta2(database_with_flanking_fasta)
    novel_allele_dict, novel_allele_name_list = parse_fasta3(corrected_alleles_raw_fasta)
    
    novel_alleles_for_retention = check_variant_pos(novel_allele_dict, novel_allele_name_list, db_without_flanking_dict, db_with_flanking_dict)
    
    with open(corrected_alleles_raw_fasta, 'w') as f:
        for line in novel_alleles_for_retention:
            f.write(line)