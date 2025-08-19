import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-dp', '--fn_depth_report',
        help = 'input read-depth calling report'
    )
    args = parser.parse_args()
    return args


def process_results(list_depth, fixed_thresh=0.5):
    total_num = 0
    novel_num = 0

    for allele_name, depth in list_depth:
        if depth > fixed_thresh:
            total_num += 1
            if 'novel' in allele_name or '_' in allele_name:
                novel_num += 1
        
        print(allele_name, depth)
        
        if depth == 0:
            return total_num, novel_num

    return total_num, novel_num



if __name__ == '__main__':
    args = parse_args()
    fn_depth_report = args.fn_depth_report

    f_n = open(fn_depth_report, 'r')
    list_depth = [] # list_depth = [(name1,depth1), (name2,depth2), ... ]
    for line in f_n:
        fields = line.split()
        allele_name = fields[0]
        depth = float(fields[1])
        list_depth.append((allele_name, depth))
        if depth == 0: # keep the first 0
            break
    f_n.close()

    total_num, novel_num = process_results(list_depth, fixed_thresh=0.5)

    print("\n========= Summary ===========")
    print("Total gAIRR-wgs alleles:", total_num)
    print("Novel gAIRR-wgs alleles:", novel_num)
    
