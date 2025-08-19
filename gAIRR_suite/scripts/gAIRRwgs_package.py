import re
def extract_mismatches(md_tag):
    matches = re.findall(r'(\d+)([A-Z]+)', md_tag)
    
    mismatch_positions = []
    current_position = 0
    
    for num, letters in matches:
        current_position += int(num)
        mismatch_positions.append(current_position + 1)  # Convert to 1-based index
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

def positions_to_check():
    # non-zero-based position
    pos_dict = {
    "TRBV21-1*01_THOL": [85, 203],
    "TRBVB*01_H7AX": [236, 379],
    "TRAV37*01_XUTF": [75, 294],
    "TRBV10-1*02": [209, 285],
    "TRBV6-8*01_XLO6": [156, 301],
    "TRBV4-2*01_IX4A": [134],
    "TRBV20/OR9-2*03": [195, 203],
    "TRAV8-4*01_NTMW": [95, 96, 204],
    "TRBV4-3*01_3IWL": [134],
    "TRBV5-8*01_RS7C": [106, 183, 287],
    "TRBV5-8*01_H3GH": [106, 183, 287],
    "TRBV20-1*01_SXO6": [78, 101, 174],
    "TRBV20-1*01_FPQG": [78, 101, 174],
    "TRBV30*02": [159, 296]
    }

    return pos_dict

def FL_region_for_check():
    gene_list = ["X04939|TRAV1-1*02|Homo","M81774|TRAV12-2*02|Homo","X04946|TRAV12-2*03|Homo","M21624|TRAV14/DV4*03|Homo","X58736|TRAV21*02|Homo",
    "X58738|TRAV35*02|Homo","Z46643|TRAV36/DV7*04|Homo","M95394|TRAV38-1*03|Homo","X58747|TRAV6*02|Homo","D13077|TRAV8-4*03|Homo","L06881|TRAV9-2*03|Homo",
    "AB305924|TRBV11-3*04|Homo","M97725|TRBV19*03|Homo","X57611|TRBV5-5*02|Homo","L36092|TRBV7-5*01|Homo","M27337|TRGV2*02|Homo","X04038|TRGV3*02|Homo",
    "Y10411|TRAV23/DV6*04|Homo","Z49060|TRAV6*03|Homo","U17051|TRBV10-1*03|Homo","L33101|TRBV10-3*03|Homo","X57722|TRBV14*02|Homo","X72719|TRBV20-1*02|Homo",
    "X57604|TRBV20-1*05|Homo","M13554|TRBV30*04|Homo","X57616|TRBV4-3*04|Homo","Y13426|TRDV2*02|Homo","AF033825|TRAJ47*02|Homo"]

    return gene_list

def get_FL_mark_dict():
    # zero-based position
    FL_mark_dict = {'new|TRAJ1*01_7NLF|Homo': [50, 111], 'new|TRAJ17*01_LF3J|Homo': [50, 112], 'new|TRAJ17*01_WFNS|Homo': [50, 112], 'new|TRAJ18*01_BZBM|Homo': [50, 115],
    'new|TRAJ19*01_VS6K|Homo': [50, 109], 'new|TRAJ3*01_FQMU|Homo': [50, 111], 'new|TRAJ36*01_C7FZ|Homo': [50, 108], 'new|TRAJ38*01_3IIW|Homo': [50, 111],
    'new|TRAJ4*01_CA4O|Homo': [50, 112], 'new|TRAJ41*01_QYDV|Homo': [50, 111], 'new|TRAJ44*01_EAA5|Homo': [50, 112], 'new|TRAJ44*01_UGG7|Homo': [50, 112],
    'new|TRAJ47*01_6LEL|Homo': [50, 106], 'new|TRAJ47*01_IZJW|Homo': [50, 106], 'new|TRAJ5*01_PVPV|Homo': [50, 109], 'new|TRAJ52*01_HKTI|Homo': [50, 118],
    'new|TRAJ52*01_JEB4|Homo': [50, 117], 'new|TRAJ53*01_ZPCD|Homo': [50, 115], 'new|TRAJ54*01_FFTX|Homo': [50, 109], 'new|TRAJ57*01_ZXI6|Homo': [50, 112],
    'new|TRAJ9*01_4QYL|Homo': [50, 110], 'new|TRBJ1-2*01_WFQ7|Homo': [50, 97], 'new|TRBJ2-3*01_FA4L|Homo': [50, 98], 'new|TRBJ2-7*01_UZLC|Homo': [50, 96],
    'new|TRDJ2*01_QSDL|Homo': [50, 103], 'new|TRDJ3*01_HT3N|Homo': [50, 108], 'new|TRDJ4*01_G4F3|Homo': [50, 97], 'new|TRGJP*01_WTZX|Homo': [50, 111],
    'new|TRGJP2*01_CUJY|Homo': [50, 109], 'new|TRGJP2*01_NXVF|Homo': [50, 109], 'new240314|TRAJ47*01_ASQO|Homo': [50, 106], 'K02545|TRBD1*01|Homo': [50, 75],
    'X02987|TRBD2*01|Homo': [50, 79], 'M14159|TRBD2*02|Homo': [50, 79], 'M23325|TRDD1*01|Homo': [50, 70], 'M22153|TRDD2*01|Homo': [50, 72],
    'M22152|TRDD3*01|Homo': [50, 76], 'X02884|TRAJ1*01|Homo': [50, 111], 'M94081|TRAJ10*01|Homo': [50, 113], 'M94081|TRAJ11*01|Homo': [50, 109],
    'X02885|TRAJ12*01|Homo': [50, 109], 'M94081|TRAJ13*01|Homo': [50, 112], 'AB258131|TRAJ13*02|Homo': [50, 112], 'M94081|TRAJ14*01|Homo': [50, 101],
    'X05775|TRAJ15*01|Homo': [50, 109], 'M94081|TRAJ15*02|Homo': [0, 59], 'M94081|TRAJ16*01|Homo': [50, 109], 'IMGT000024|TRAJ16*02|Homo': [50, 109],
    'X05773|TRAJ17*01|Homo': [50, 112], 'M94081|TRAJ18*01|Homo': [50, 115], 'M94081|TRAJ19*01|Homo': [50, 109], 'X02884|TRAJ2*01|Homo': [50, 115],
    'M94081|TRAJ20*01|Homo': [50, 106], 'M94081|TRAJ21*01|Homo': [50, 104], 'X02886|TRAJ22*01|Homo': [50, 112], 'M94081|TRAJ23*01|Homo': [50, 112],
    'X58763|TRAJ23*02|Homo': [0, 62], 'X02887|TRAJ24*01|Homo': [0, 62], 'M94081|TRAJ24*02|Homo': [50, 112], 'IMGT000024|TRAJ24*03|Homo': [50, 112],
    'X02888|TRAJ25*01|Homo': [50, 109], 'M94081|TRAJ26*01|Homo': [50, 109], 'M94081|TRAJ27*01|Homo': [50, 108], 'M94081|TRAJ28*01|Homo': [50, 115],
    'M94081|TRAJ29*01|Homo': [50, 109], 'X02884|TRAJ3*01|Homo': [50, 111], 'M94081|TRAJ30*01|Homo': [50, 106], 'M94081|TRAJ31*01|Homo': [50, 106],
    'M94081|TRAJ32*01|Homo': [50, 115], 'IMGT000024|TRAJ32*02|Homo': [50, 115], 'M94081|TRAJ33*01|Homo': [50, 106], 'M35622|TRAJ34*01|Homo': [50, 107],
    'M94081|TRAJ35*01|Homo': [50, 108], 'M94081|TRAJ36*01|Homo': [50, 108], 'M94081|TRAJ37*01|Homo': [50, 111], 'IMGT000024|TRAJ37*02|Homo': [50, 111],
    'M94081|TRAJ38*01|Homo': [50, 111], 'M94081|TRAJ39*01|Homo': [50, 112], 'M94081|TRAJ4*01|Homo': [50, 112], 'M35620|TRAJ40*01|Homo': [50, 110],
    'M94081|TRAJ41*01|Homo': [50, 111], 'M94081|TRAJ42*01|Homo': [50, 115], 'M94081|TRAJ43*01|Homo': [50, 103], 'M35619|TRAJ44*01|Homo': [50, 112],
    'M94081|TRAJ45*01|Homo': [50, 115], 'M94081|TRAJ46*01|Homo': [50, 112], 'M94081|TRAJ47*01|Homo': [50, 106], 'AF033825|TRAJ47*02|Homo': [50, 105],
    'M94081|TRAJ48*01|Homo': [50, 112], 'M94081|TRAJ49*01|Homo': [50, 105], 'M94081|TRAJ5*01|Homo': [50, 109], 'M94081|TRAJ50*01|Homo': [50, 109],
    'M94081|TRAJ51*01|Homo': [50, 112], 'M94081|TRAJ52*01|Homo': [50, 118], 'M94081|TRAJ53*01|Homo': [50, 115], 'M94081|TRAJ54*01|Homo': [50, 109],
    'M94081|TRAJ55*01|Homo': [50, 106], 'M94081|TRAJ56*01|Homo': [50, 111], 'M94081|TRAJ57*01|Homo': [50, 112], 'M94081|TRAJ58*01|Homo': [50, 112],
    'M94081|TRAJ59*01|Homo': [50, 109], 'M16747|TRAJ6*01|Homo': [50, 111], 'M94081|TRAJ60*01|Homo': [50, 106], 'M94081|TRAJ61*01|Homo': [50, 109],
    'M94081|TRAJ7*01|Homo': [50, 108], 'M94081|TRAJ8*01|Homo': [50, 109], 'IMGT000024|TRAJ8*02|Homo': [0, 59], 'M94081|TRAJ9*01|Homo': [50, 110],
    'K02545|TRBJ1-1*01|Homo': [50, 97], 'K02545|TRBJ1-2*01|Homo': [50, 97], 'M14158|TRBJ1-3*01|Homo': [50, 99], 'M14158|TRBJ1-4*01|Homo': [50, 100],
    'M14158|TRBJ1-5*01|Homo': [50, 99], 'M14158|TRBJ1-6*01|Homo': [50, 102], 'L36092|TRBJ1-6*02|Homo': [50, 102], 'X02987|TRBJ2-1*01|Homo': [50, 99],
    'X02987|TRBJ2-2*01|Homo': [50, 100], 'X02987|TRBJ2-2P*01|Homo': [50, 95], 'X02987|TRBJ2-3*01|Homo': [50, 98], 'X02987|TRBJ2-4*01|Homo': [50, 99],
    'X02987|TRBJ2-5*01|Homo': [50, 97], 'X02987|TRBJ2-6*01|Homo': [50, 102], 'M14159|TRBJ2-7*01|Homo': [50, 96], 'X02987|TRBJ2-7*02|Homo': [50, 96],
    'M20289|TRDJ1*01|Homo': [50, 100], 'L36386|TRDJ2*01|Homo': [50, 103], 'M21508|TRDJ3*01|Homo': [50, 108], 'AJ249814|TRDJ4*01|Homo': [50, 97],
    'M12960|TRGJ1*01|Homo': [50, 99], 'IMGT000011|TRGJ1*02|Homo': [50, 99], 'M12961|TRGJ2*01|Homo': [50, 99], 'M12950|TRGJP*01|Homo': [50, 111],
    'X08084|TRGJP1*01|Homo': [50, 109], 'M16016|TRGJP2*01|Homo': [50, 109], 'AE000658|TRAV1-1*01|Homo': [50, 324], 'X04939|TRAV1-1*02|Homo': [50, 318],
    'AE000658|TRAV1-2*01|Homo': [50, 324], 'U32544|TRAV1-2*02|Homo': [0, 175], 'IMGT000024|TRAV1-2*03|Homo': [50, 324], 'AE000659|TRAV10*01|Homo': [50, 329],
    'IMGT000024|TRAV10*02|Homo': [50, 329], 'AE000659|TRAV11*01|Homo': [50, 326], 'AE000659|TRAV11-1*01|Homo': [50, 329], 'AE000659|TRAV12-1*01|Homo': [50, 323],
    'M17657|TRAV12-1*02|Homo': [0, 273], 'AE000659|TRAV12-2*01|Homo': [50, 326], 'M81774|TRAV12-2*02|Homo': [50, 322], 'X04946|TRAV12-2*03|Homo': [50, 301],
    'AE000659|TRAV12-3*01|Homo': [50, 326], 'M17656|TRAV12-3*02|Homo': [0, 276], 'AE000659|TRAV13-1*01|Homo': [50, 329], 'X04954|TRAV13-1*02|Homo': [50, 329],
    'L11162|TRAV13-1*03|Homo': [0, 233], 'AE000659|TRAV13-2*01|Homo': [50, 329], 'M17658|TRAV13-2*02|Homo': [0, 276], 'AE000659|TRAV14-1*01|Homo': [50, 333],
    'IMGT000024|TRAV14-1*02|Homo': [50, 333], 'M21626|TRAV14/DV4*01|Homo': [50, 339], 'AE000659|TRAV14/DV4*02|Homo': [50, 339], 'M21624|TRAV14/DV4*03|Homo': [50, 331],
    'L09758|TRAV14/DV4*04|Homo': [0, 270], 'AE000659|TRAV15*01|Homo': [50, 332], 'AE000659|TRAV16*01|Homo': [50, 321], 'AE000660|TRAV17*01|Homo': [50, 326],
    'AE000660|TRAV18*01|Homo': [50, 327], 'AE000660|TRAV19*01|Homo': [50, 339], 'AE000658|TRAV2*01|Homo': [50, 312], 'M17659|TRAV2*02|Homo': [0, 263],
    'AE000660|TRAV20*01|Homo': [50, 323], 'IMGT000024|TRAV20*02|Homo': [50, 323], 'S60789|TRAV20*03|Homo': [0, 263], 'X70305|TRAV20*04|Homo': [0, 266],
    'AE000660|TRAV21*01|Homo': [50, 328], 'X58736|TRAV21*02|Homo': [50, 322], 'AE000660|TRAV22*01|Homo': [50, 320], 'AE000660|TRAV23/DV6*01|Homo': [50, 329],
    'M17660|TRAV23/DV6*02|Homo': [0, 279], 'M97704|TRAV23/DV6*03|Homo': [50, 329], 'Y10411|TRAV23/DV6*04|Homo': [50, 316], 'IMGT000024|TRAV23/DV6*05|Homo': [0, 279],
    'AE000660|TRAV24*01|Homo': [50, 326], 'new|TRAV24*01_TWOO|Homo': [50, 326],'new|TRAV24*01_THRE|Homo': [50, 326],'new|TRAV24*01_FOUR|Homo': [50, 326],'M17661|TRAV24*02|Homo': [0, 276], 'AE000660|TRAV25*01|Homo': [50, 322], 'AE000660|TRAV26-1*01|Homo': [50, 326],
    'IMGT000024|TRAV26-1*02|Homo': [50, 326], 'L06886|TRAV26-1*03|Homo': [0, 266], 'AE000660|TRAV26-2*01|Homo': [50, 325], 'L11160|TRAV26-2*02|Homo': [0, 216],
    'AE000660|TRAV27*01|Homo': [50, 323], 'X04957|TRAV27*02|Homo': [0, 271], 'IMGT000024|TRAV27*03|Homo': [50, 323], 'AE000660|TRAV28*01|Homo': [50, 318],
    'IMGT000024|TRAV28*02|Homo': [50, 318], 'AE000660|TRAV29/DV5*01|Homo': [50, 329], 'S81645|TRAV29/DV5*02|Homo': [0, 278], 'AJ245565|TRAV29/DV5*03|Homo': [50, 328],
    'IMGT000024|TRAV29/DV5*04|Homo': [50, 329], 'AE000658|TRAV3*01|Homo': [50, 335], 'M27377|TRAV3*02|Homo': [0, 283], 'AE000660|TRAV30*01|Homo': [50, 323],
    'X58768|TRAV30*02|Homo': [0, 266], 'L06883|TRAV30*03|Homo': [0, 266], 'U32537|TRAV30*04|Homo': [0, 229], 'IMGT000024|TRAV30*05|Homo': [50, 323],
    'AE000660|TRAV31*01|Homo': [50, 335], 'IMGT000024|TRAV31*02|Homo': [50, 335], 'AE000660|TRAV32*01|Homo': [50, 325], 'AE000660|TRAV33*01|Homo': [50, 327],
    'Z46653|TRAV33*02|Homo': [0, 268], 'AE000660|TRAV34*01|Homo': [50, 326], 'AE000660|TRAV35*01|Homo': [50, 325], 'X58738|TRAV35*02|Homo': [50, 319],
    'IMGT000024|TRAV35*03|Homo': [50, 325], 'AE000660|TRAV36/DV7*01|Homo': [50, 326], 'X61070|TRAV36/DV7*02|Homo': [0, 273], 'X58767|TRAV36/DV7*03|Homo': [0, 269],
    'Z46643|TRAV36/DV7*04|Homo': [50, 320], 'IMGT000024|TRAV36/DV7*05|Homo': [50, 326], 'AE000661|TRAV37*01|Homo': [50, 328], 'AE000661|TRAV38-1*01|Homo': [50, 339],
    'M64355|TRAV38-1*02|Homo': [0, 279], 'M95394|TRAV38-1*03|Homo': [50, 332], 'L06880|TRAV38-1*04|Homo': [0, 278], 'AE000661|TRAV38-2/DV8*01|Homo': [50, 338],
    'AE000661|TRAV39*01|Homo': [50, 326], 'AE000658|TRAV4*01|Homo': [50, 326], 'X73521|TRAV40*01|Homo': [50, 312], 'AE000661|TRAV41*01|Homo': [50, 319],
    'AE000661|TRAV46*01|Homo': [50, 316], 'AE000659|TRAV5*01|Homo': [50, 326], 'AE000659|TRAV6*01|Homo': [50, 329], 'X58747|TRAV6*02|Homo': [50, 322],
    'Z49060|TRAV6*03|Homo': [50, 298], 'Y10409|TRAV6*04|Homo': [0, 248], 'Y10410|TRAV6*05|Homo': [0, 248], 'U32542|TRAV6*06|Homo': [0, 216],
    'IMGT000024|TRAV6*07|Homo': [50, 329], 'AE000659|TRAV7*01|Homo': [50, 323], 'AE000659|TRAV8-1*01|Homo': [50, 333], 'U32520|TRAV8-1*02|Homo': [0, 249],
    'AE000659|TRAV8-2*01|Homo': [50, 333], 'M17650|TRAV8-2*02|Homo': [0, 279], 'IMGT000024|TRAV8-2*03|Homo': [50, 333], 'AE000659|TRAV8-3*01|Homo': [50, 333],
    'M35617|TRAV8-3*02|Homo': [0, 281], 'L06885|TRAV8-3*03|Homo': [0, 275], 'AE000659|TRAV8-4*01|Homo': [50, 333], 'M12423|TRAV8-4*02|Homo': [0, 283],
    'D13077|TRAV8-4*03|Homo': [50, 325], 'M12959|TRAV8-4*04|Homo': [0, 283], 'X63455|TRAV8-4*05|Homo': [0, 283], 'K02777|TRAV8-4*06|Homo': [0, 217],
    'M17665|TRAV8-4*07|Homo': [0, 196], 'AE000659|TRAV8-5*01|Homo': [50, 133], 'X02850|TRAV8-6*01|Homo': [50, 333], 'AE000659|TRAV8-6*02|Homo': [50, 333],
    'AE000660|TRAV8-6-1*01|Homo': [50, 335], 'AE000660|TRAV8-7*01|Homo': [50, 333], 'U76492|TRAV8-7*02|Homo': [50, 333], 'AE000659|TRAV9-1*01|Homo': [50, 330],
    'AE000659|TRAV9-2*01|Homo': [50, 330], 'IMGT000024|TRAV9-2*02|Homo': [50, 330], 'L06881|TRAV9-2*03|Homo': [50, 322], 'L06882|TRAV9-2*04|Homo': [0, 272],
    'AE000659|TRAVA*01|Homo': [50, 334], 'IMGT000024|TRAVA*02|Homo': [50, 334], 'AE000659|TRAVB*01|Homo': [50, 317], 'IMGT000024|TRAVB*02|Homo': [50, 319],
    'AE000660|TRAVC*01|Homo': [50, 337], 'L36092|TRBV1*01|Homo': [50, 333], 'L36092|TRBV10-1*01|Homo': [50, 336], 'AF009660|TRBV10-1*02|Homo': [50, 336],
    'U17051|TRBV10-1*03|Homo': [50, 264], 'L36092|TRBV10-2*01|Homo': [50, 336], 'IMGT000021|TRBV10-2*02|Homo': [50, 336], 'U03115|TRBV10-3*01|Homo': [50, 336],
    'U17047|TRBV10-3*02|Homo': [50, 336], 'L33101|TRBV10-3*03|Homo': [50, 322], 'L33102|TRBV10-3*04|Homo': [0, 272], 'M33233|TRBV11-1*01|Homo': [50, 339],
    'L36092|TRBV11-2*01|Homo': [50, 339], 'M33235|TRBV11-2*02|Homo': [0, 284], 'IMGT000021|TRBV11-2*03|Homo': [50, 339], 'U03115|TRBV11-3*01|Homo': [50, 339],
    'X58797|TRBV11-3*02|Homo': [0, 284], 'M62377|TRBV11-3*03|Homo': [0, 268], 'AB305924|TRBV11-3*04|Homo': [50, 338], 'X07224|TRBV12-1*01|Homo': [50, 339],
    'X06936|TRBV12-2*01|Homo': [50, 339], 'X07192|TRBV12-3*01|Homo': [50, 339], 'K02546|TRBV12-4*01|Homo': [50, 339], 'M14264|TRBV12-4*02|Homo': [0, 287],
    'X07223|TRBV12-5*01|Homo': [50, 339], 'U03115|TRBV13*01|Homo': [50, 336], 'M62378|TRBV13*02|Homo': [0, 281], 'X06154|TRBV14*01|Homo': [50, 339],
    'X57722|TRBV14*02|Homo': [50, 334], 'U03115|TRBV15*01|Homo': [50, 336], 'IMGT000021|TRBV15*02|Homo': [50, 336], 'M62376|TRBV15*03|Homo': [0, 281],
    'L26231|TRBV16*01|Homo': [50, 339], 'U03115|TRBV16*02|Homo': [50, 339], 'L26054|TRBV16*03|Homo': [0, 284], 'U03115|TRBV17*01|Homo': [50, 336],
    'IMGT000021|TRBV17*02|Homo': [50, 336], 'L36092|TRBV18*01|Homo': [50, 339], 'L36092|TRBV19*01|Homo': [50, 336], 'U48259|TRBV19*02|Homo': [0, 286],
    'M97725|TRBV19*03|Homo': [50, 330], 'L36092|TRBV2*01|Homo': [50, 339], 'M62379|TRBV2*02|Homo': [50, 334], 'M64351|TRBV2*03|Homo': [0, 287],
    'M11955|TRBV20-1*01|Homo': [50, 342], 'X72719|TRBV20-1*02|Homo': [50, 337], 'M11954|TRBV20-1*03|Homo': [0, 287], 'M14263|TRBV20-1*04|Homo': [0, 290],
    'X57604|TRBV20-1*05|Homo': [50, 340], 'D13088|TRBV20-1*06|Homo': [0, 287], 'X74852|TRBV20-1*07|Homo': [0, 290], 'L05149|TRBV20/OR9-2*01|Homo': [50, 342],
    'X72718|TRBV20/OR9-2*02|Homo': [0, 287], 'AF029308|TRBV20/OR9-2*03|Homo': [50, 342], 'L36092|TRBV21-1*01|Homo': [50, 339], 'IMGT000021|TRBV21-1*02|Homo': [50, 339],
    'AF029308|TRBV21/OR9-2*01|Homo': [50, 339], 'L36092|TRBV22-1*01|Homo': [50, 337], 'AF029308|TRBV22/OR9-2*01|Homo': [50, 336], 'L36092|TRBV23-1*01|Homo': [50, 339],
    'AF029308|TRBV23/OR9-2*01|Homo': [50, 339], 'L27615|TRBV23/OR9-2*02|Homo': [0, 230], 'M11951|TRBV24-1*01|Homo': [50, 337], 'IMGT000021|TRBV24-1*02|Homo': [50, 337],
    'L05153|TRBV24/OR9-2*01|Homo': [50, 337], 'L27613|TRBV24/OR9-2*02|Homo': [0, 228], 'AF029308|TRBV24/OR9-2*03|Homo': [50, 342], 'L36092|TRBV25-1*01|Homo': [50, 336],
    'L05152|TRBV25/OR9-2*01|Homo': [50, 330], 'L27611|TRBV25/OR9-2*02|Homo': [0, 245], 'L36092|TRBV26*01|Homo': [50, 336], 'AF029308|TRBV26/OR9-2*01|Homo': [50, 336],
    'AL356489|TRBV26/OR9-2*02|Homo': [50, 336], 'L36092|TRBV27*01|Homo': [50, 336], 'U08314|TRBV28*01|Homo': [50, 336], 'L36092|TRBV29-1*01|Homo': [50, 339],
    'M13847|TRBV29-1*02|Homo': [0, 287], 'X04926|TRBV29-1*03|Homo': [0, 230], 'L05150|TRBV29/OR9-2*01|Homo': [0, 289], 'AF029308|TRBV29/OR9-2*02|Homo': [50, 339],
    'U07977|TRBV3-1*01|Homo': [50, 336], 'L06889|TRBV3-1*02|Homo': [0, 278], 'L36092|TRBV3-2*01|Homo': [50, 336], 'U07978|TRBV3-2*02|Homo': [50, 336],
    'M33240|TRBV3-2*03|Homo': [0, 285], 'L36092|TRBV30*01|Homo': [50, 333], 'Z13967|TRBV30*02|Homo': [50, 333], 'IMGT000027|TRBV30*03|Homo': [0, 283],
    'M13554|TRBV30*04|Homo': [50, 325], 'L06893|TRBV30*05|Homo': [0, 281], 'U07977|TRBV4-1*01|Homo': [50, 336], 'M13855|TRBV4-1*02|Homo': [0, 258],
    'U07975|TRBV4-2*01|Homo': [50, 336], 'X58811|TRBV4-2*02|Homo': [0, 281], 'U07978|TRBV4-3*01|Homo': [50, 336], 'X58812|TRBV4-3*02|Homo': [0, 281],
    'L06888|TRBV4-3*03|Homo': [0, 281], 'X57616|TRBV4-3*04|Homo': [50, 280], 'L36092|TRBV5-1*01|Homo': [50, 335], 'M14271|TRBV5-1*02|Homo': [0, 278],
    'L36092|TRBV5-2*01|Homo': [50, 326], 'X61439|TRBV5-3*01|Homo': [50, 335], 'AF009660|TRBV5-3*02|Homo': [50, 335], 'L36092|TRBV5-4*01|Homo': [50, 335],
    'X57615|TRBV5-4*02|Homo': [0, 281], 'S50547|TRBV5-4*03|Homo': [0, 233], 'X58804|TRBV5-4*04|Homo': [0, 191], 'L36092|TRBV5-5*01|Homo': [50, 335],
    'X57611|TRBV5-5*02|Homo': [50, 331], 'X58801|TRBV5-5*03|Homo': [0, 281], 'L36092|TRBV5-6*01|Homo': [50, 335], 'L36092|TRBV5-7*01|Homo': [50, 335],
    'L36092|TRBV5-8*01|Homo': [50, 335], 'X58803|TRBV5-8*02|Homo': [0, 237], 'X61446|TRBV6-1*01|Homo': [50, 336], 'X61445|TRBV6-2*01|Homo': [50, 336],
    'U07978|TRBV6-3*01|Homo': [50, 336], 'X61653|TRBV6-4*01|Homo': [50, 336], 'AF009660|TRBV6-4*02|Homo': [50, 336], 'L36092|TRBV6-5*01|Homo': [50, 336],
    'L36092|TRBV6-6*01|Homo': [50, 336], 'AF009662|TRBV6-6*02|Homo': [50, 336], 'X58815|TRBV6-6*03|Homo': [0, 281], 'X74848|TRBV6-6*04|Homo': [0, 284],
    'L06892|TRBV6-6*05|Homo': [0, 281], 'L36092|TRBV6-7*01|Homo': [50, 336], 'L36092|TRBV6-8*01|Homo': [50, 333], 'X61447|TRBV6-9*01|Homo': [50, 336],
    'X61444|TRBV7-1*01|Homo': [50, 339], 'X61442|TRBV7-2*01|Homo': [50, 339], 'L36190|TRBV7-2*02|Homo': [50, 339], 'U07975|TRBV7-2*03|Homo': [0, 289],
    'M27387|TRBV7-2*04|Homo': [0, 288], 'X61440|TRBV7-3*01|Homo': [50, 339], 'M97943|TRBV7-3*02|Homo': [50, 339], 'AF009660|TRBV7-3*03|Homo': [50, 339],
    'X74843|TRBV7-3*04|Homo': [0, 286], 'M13550|TRBV7-3*05|Homo': [0, 230], 'L36092|TRBV7-4*01|Homo': [50, 339], 'L13762|TRBV7-4*02|Homo': [0, 259],
    'L36092|TRBV7-5*01|Homo': [50, 337], 'AF009663|TRBV7-5*02|Homo': [50, 338], 'L36092|TRBV7-6*01|Homo': [50, 339], 'X58806|TRBV7-6*02|Homo': [0, 284],
    'L36092|TRBV7-7*01|Homo': [50, 339], 'X57607|TRBV7-7*02|Homo': [0, 284], 'M11953|TRBV7-8*01|Homo': [50, 339], 'X61441|TRBV7-8*02|Homo': [50, 339],
    'M27384|TRBV7-8*03|Homo': [0, 287], 'L36092|TRBV7-9*01|Homo': [50, 339], 'M15564|TRBV7-9*02|Homo': [0, 289], 'AF009663|TRBV7-9*03|Homo': [50, 339],
    'M14261|TRBV7-9*04|Homo': [0, 284], 'M27385|TRBV7-9*05|Homo': [0, 287], 'X74844|TRBV7-9*06|Homo': [0, 287], 'L14854|TRBV7-9*07|Homo': [0, 202],
    'L36092|TRBV8-1*01|Homo': [50, 328], 'IMGT000021|TRBV8-1*02|Homo': [50, 328], 'L36092|TRBV8-2*01|Homo': [50, 322], 'IMGT000021|TRBV8-2*02|Homo': [50, 322],
    'L36092|TRBV9*01|Homo': [50, 335], 'AF009660|TRBV9*02|Homo': [50, 335], 'M27380|TRBV9*03|Homo': [0, 281], 'L36092|TRBVA*01|Homo': [50, 311],
    'IMGT000021|TRBVA*02|Homo': [50, 311], 'AF029308|TRBVA/OR9-2*01|Homo': [50, 311], 'L36092|TRBVB*01|Homo': [50, 401], 'IMGT000021|TRBVB*02|Homo': [50, 401],
    'L36092|TRBVC*01|Homo': [50, 138], 'M22198|TRDV1*01|Homo': [50, 336], 'X15207|TRDV2*01|Homo': [50, 337], 'Y13426|TRDV2*02|Homo': [50, 332],
    'AE000661|TRDV2*03|Homo': [50, 337], 'M23326|TRDV3*01|Homo': [50, 339], 'X15261|TRDV3*02|Homo': [50, 339], 'M12949|TRGV1*01|Homo': [50, 346],
    'X07206|TRGV10*01|Homo': [0, 310], 'X74798|TRGV10*02|Homo': [50, 358], 'Y11227|TRGV11*01|Homo': [50, 358], 'AF159056|TRGV11*02|Homo': [0, 308],
    'M13429|TRGV2*01|Homo': [50, 349], 'M27337|TRGV2*02|Homo': [50, 346], 'IMGT000011|TRGV2*03|Homo': [50, 349], 'M13430|TRGV3*01|Homo': [50, 349],
    'X04038|TRGV3*02|Homo': [50, 348], 'X15272|TRGV4*01|Homo': [50, 349], 'X13354|TRGV4*02|Homo': [50, 349], 'X13355|TRGV5*01|Homo': [50, 349],
    'M13431|TRGV5P*01|Homo': [0, 299], 'AF159056|TRGV5P*02|Homo': [50, 349], 'M13432|TRGV6*01|Homo': [50, 348], 'AF159056|TRGV6*02|Homo': [50, 348],
    'M13433|TRGV7*01|Homo': [50, 348], 'M13434|TRGV8*01|Homo': [50, 349], 'X07205|TRGV9*01|Homo': [50, 355], 'X15274|TRGV9*02|Homo': [0, 305],
    'X07208|TRGVA*01|Homo': [50, 334], 'X07209|TRGVB*01|Homo': [50, 361], 'IMGT000011|TRGVB*02|Homo': [50, 361], 'new|TRAV10*01_DOMV|Homo': [50, 329],
    'new|TRAV10*01_IM77|Homo': [50, 329], 'new|TRAV10*01_S23D|Homo': [50, 329], 'new|TRAV11*01_B7E5|Homo': [50, 326], 'new|TRAV11*01_P4P2|Homo': [50, 326],
    'new|TRAV11*01_PCE4|Homo': [50, 326], 'new|TRAV1-1*01_NHW2|Homo': [50, 324], 'new|TRAV1-1*01_YVAF|Homo': [50, 324], 'new|TRAV11-1*01_CWDU|Homo': [50, 329],
    'new|TRAV1-2*03_D5LG|Homo': [50, 323], 'new|TRAV12-1*01_PA6M|Homo': [50, 323], 'new|TRAV12-1*01_UE5Z|Homo': [50, 323], 'new|TRAV12-1*01_XNKU|Homo': [50, 323],
    'new|TRAV12-2*01_7QKH|Homo': [50, 326], 'new|TRAV12-3*01_XDI7|Homo': [50, 326], 'new|TRAV13-1*02_IFSB|Homo': [50, 329], 'new|TRAV13-2*01_3UGI|Homo': [50, 329],
    'new|TRAV13-2*01_CCFZ|Homo': [50, 329], 'new|TRAV13-2*01_HHMW|Homo': [50, 329], 'new|TRAV14/DV4*01_4XZE|Homo': [50, 339], 'new|TRAV14/DV4*02_COXU|Homo': [50, 339],
    'new|TRAV14/DV4*02_TAKJ|Homo': [50, 339], 'new|TRAV14-1*01_46TF|Homo': [50, 333], 'new|TRAV14-1*01_R26U|Homo': [50, 333], 'new|TRAV15*01_26JH|Homo': [50, 332],
    'new|TRAV15*01_HNSX|Homo': [50, 332], 'new|TRAV15*01_JNCB|Homo': [50, 332], 'new|TRAV15*01_PKAX|Homo': [50, 332], 'new|TRAV15*01_UCH6|Homo': [50, 332],
    'new|TRAV18*01_HOW3|Homo': [50, 327], 'new|TRAV19*01_MTHX|Homo': [50, 339], 'new|TRAV19*01_VFBW|Homo': [50, 339], 'new|TRAV19*01_XX6B|Homo': [50, 339],
    'new|TRAV19*01_YWUQ|Homo': [50, 339], 'new|TRAV2*01_GJQM|Homo': [50, 312], 'new|TRAV23/DV6*01_BQK2|Homo': [50, 329], 'new|TRAV24*01_JERO|Homo': [50, 326],
    'new|TRAV25*01_JW7E|Homo': [50, 322], 'new|TRAV25*01_QZWI|Homo': [50, 321], 'new|TRAV25*01_ZHIK|Homo': [50, 322], 'new|TRAV26-1*01_BZTF|Homo': [50, 327],
    'new|TRAV26-2*01_5R7N|Homo': [50, 325], 'new|TRAV27*01_MPWQ|Homo': [50, 323], 'new|TRAV27*03_JHU5|Homo': [50, 323], 'new|TRAV30*05_MCUK|Homo': [50, 323],
    'new|TRAV30*05_XW5Q|Homo': [50, 323], 'new|TRAV31*01_3EPX|Homo': [50, 335], 'new|TRAV31*01_3M4G|Homo': [50, 335], 'new|TRAV31*01_QVZD|Homo': [50, 335],
    'new|TRAV31*01_WEDB|Homo': [50, 335], 'new|TRAV31*01_3EXU|Homo': [50, 335], 'new|TRAV32*01_2A4Y|Homo': [50, 325], 'new|TRAV32*01_ALWP|Homo': [50, 325],
    'new|TRAV32*01_C5YN|Homo': [50, 325], 'new|TRAV32*01_X3HG|Homo': [50, 325], 'new|TRAV33*01_MBZQ|Homo': [50, 327], 'new|TRAV34*01_5GW5|Homo': [50, 326],
    'new|TRAV37*01_Q43H|Homo': [50, 328], 'new|TRAV37*01_XUTF|Homo': [50, 328], 'new|TRAV4*01_OEZV|Homo': [50, 326], 'new|TRAV40*01_O36H|Homo': [50, 312],
    'new|TRAV40*01_POMY|Homo': [50, 312], 'new|TRAV41*01_ENO4|Homo': [50, 318], 'new|TRAV41*01_PUIX|Homo': [50, 319], 'new|TRAV46*01_DDKN|Homo': [50, 316],
    'new|TRAV7*01_PYBO|Homo': [50, 323], 'new|TRAV8-1*01_JVON|Homo': [50, 333], 'new|TRAV8-1*01_P2AT|Homo': [50, 333], 'new|TRAV8-1*01_PEFX|Homo': [50, 333],
    'new|TRAV8-2*03_2GPR|Homo': [50, 333], 'new|TRAV8-3*01_RZJS|Homo': [50, 333], 'new|TRAV8-3*01_YUIS|Homo': [50, 333], 'new|TRAV8-4*01_GWC4|Homo': [50, 333],
    'new|TRAV8-4*01_NTMW|Homo': [50, 333], 'new|TRAV8-4*01_ODFF|Homo': [50, 333], 'new|TRAV8-4*01_Q6HW|Homo': [50, 333], 'new|TRAV8-4*01_A6EG|Homo': [50, 333],
    'new|TRAV8-5*01_2IYY|Homo': [50, 133], 'new|TRAV8-6*01_X6OM|Homo': [50, 333], 'new|TRAV8-6-1*01_6OTT|Homo': [50, 335], 'new|TRAV8-6-1*01_76AO|Homo': [50, 335],
    'new|TRAV8-6-1*01_JXVH|Homo': [50, 335], 'new|TRAV8-6-1*01_ODNT|Homo': [50, 335], 'new|TRAV9-1*01_TOQ5|Homo': [50, 330], 'new|TRAV9-1*01_ZCB6|Homo': [50, 330],
    'new|TRAV9-2*02_27EV|Homo': [50, 330], 'new|TRAVA*01_GFYT|Homo': [50, 334], 'new|TRAVA*02_4BST|Homo': [50, 334], 'new|TRAVA*02_5FP6|Homo': [50, 335],
    'new|TRAVA*02_FSXW|Homo': [50, 334], 'new|TRAVA*02_HE43|Homo': [50, 334], 'new|TRAVC*01_JUHJ|Homo': [50, 337], 'new|TRBV1*01_ES47|Homo': [50, 333],
    'new|TRBV10-2*02_BZMU|Homo': [50, 336], 'new|TRBV11-1*01_B6R4|Homo': [50, 340], 'new|TRBV11-1*01_ZEZ2|Homo': [50, 339], 'new|TRBV11-2*01_SBUO|Homo': [50, 325],
    'new|TRBV11-2*01_U4GQ|Homo': [50, 339], 'new|TRBV12-1*01_K7L2|Homo': [50, 339], 'new|TRBV12-1*01_Q6AZ|Homo': [50, 339], 'new|TRBV12-1*01_QZMG|Homo': [50, 339],
    'new|TRBV12-2*01_6S66|Homo': [50, 339], 'new|TRBV12-2*01_T745|Homo': [50, 339], 'new|TRBV12-3*01_7SWL|Homo': [50, 339], 'new|TRBV12-3*01_HYTG|Homo': [50, 339],
    'new|TRBV12-4*01_TUFA|Homo': [50, 339], 'new|TRBV12-4*01_UNNP|Homo': [50, 339], 'new|TRBV12-4*01_WZEK|Homo': [50, 339], 'new|TRBV12-4*01_YECY|Homo': [50, 339],
    'new|TRBV12-5*01_C3AM|Homo': [50, 339], 'new|TRBV12-5*01_GMUN|Homo': [50, 339], 'new|TRBV13*01_ATEA|Homo': [50, 335], 'new|TRBV13*01_FUIZ|Homo': [50, 336],
    'new|TRBV13*01_UN4K|Homo': [50, 336], 'new|TRBV14*01_4BXP|Homo': [50, 339], 'new|TRBV15*02_AO6X|Homo': [50, 336], 'new|TRBV15*02_OGKD|Homo': [50, 336],
    'new|TRBV16*02_7LFA|Homo': [50, 339], 'new|TRBV17*01_XV5U|Homo': [50, 336], 'new|TRBV18*01_HWMM|Homo': [50, 339], 'new|TRBV18*01_XUQB|Homo': [50, 339],
    'new|TRBV19*01_3P3P|Homo': [50, 336], 'new|TRBV20/OR9-2*01_POQ5|Homo': [50, 342], 'new|TRBV20-1*01_FPQG|Homo': [50, 342], 'new|TRBV21-1*01_IKOP|Homo': [50, 339],
    'new|TRBV21-1*01_THOL|Homo': [50, 339], 'new|TRBV21-1*01_HF7P|Homo': [50, 339], 'new|TRBV22/OR9-2*01_JGK6|Homo': [50, 336], 'new|TRBV22-1*01_32CZ|Homo': [50, 337],
    'new|TRBV22-1*01_NKOP|Homo': [50, 337], 'new|TRBV22-1*01_T2R4|Homo': [50, 337], 'new|TRBV24/OR9-2*03_46YL|Homo': [50, 342], 'new|TRBV24/OR9-2*03_CMR7|Homo': [50, 342],
    'new|TRBV24-1*01_2X3K|Homo': [50, 337], 'new|TRBV25/OR9-2*01_KVE4|Homo': [50, 330], 'new|TRBV25/OR9-2*01_RMXL|Homo': [50, 330], 'new|TRBV25-1*01_HE6H|Homo': [50, 336],
    'new|TRBV25-1*01_Y3JN|Homo': [50, 336], 'new|TRBV26*01_FSE6|Homo': [50, 336], 'new|TRBV26*01_PF3L|Homo': [50, 336], 'new|TRBV26*01_ZDDG|Homo': [50, 336],
    'new|TRBV26/OR9-2*02_2XIP|Homo': [50, 336], 'new|TRBV26/OR9-2*02_R5UD|Homo': [50, 336], 'new|TRBV27*01_5P4Q|Homo': [50, 336], 'new|TRBV27*01_CUDO|Homo': [50, 336],
    'new|TRBV27*01_OF43|Homo': [50, 336], 'new|TRBV28*01_6SPE|Homo': [50, 336], 'new|TRBV29-1*01_7SQU|Homo': [50, 339], 'new|TRBV29-1*01_CZKN|Homo': [50, 339],
    'new|TRBV30*03_3Q4G|Homo': [0, 283], 'new|TRBV30*03_RTTI|Homo': [0, 283], 'new|TRBV3-1*01_HLFN|Homo': [50, 336], 'new|TRBV3-2*01_LTA6|Homo': [50, 336],
    'new|TRBV4-1*01_CKZP|Homo': [50, 336], 'new|TRBV4-1*01_KIFG|Homo': [50, 336], 'new|TRBV4-2*01_IX4A|Homo': [50, 336], 'new|TRBV5-1*01_NCWC|Homo': [50, 335],
    'new|TRBV5-1*01_ZAVT|Homo': [50, 334], 'new|TRBV5-2*01_764S|Homo': [50, 326], 'new|TRBV5-2*01_QOXS|Homo': [50, 326], 'new|TRBV5-2*01_WT7U|Homo': [50, 326],
    'new|TRBV5-3*01_KVLZ|Homo': [50, 335], 'new|TRBV5-3*01_PUP4|Homo': [50, 335], 'new|TRBV5-4*01_YASX|Homo': [50, 335], 'new|TRBV5-5*01_A74Z|Homo': [50, 335],
    'new|TRBV5-5*01_CAMA|Homo': [50, 335], 'new|TRBV5-6*01_DPX7|Homo': [50, 335], 'new|TRBV5-6*01_FCU3|Homo': [50, 335], 'new|TRBV5-6*01_LGQ4|Homo': [50, 335],
    'new|TRBV5-6*01_T7HS|Homo': [50, 335], 'new|TRBV5-7*01_7OVA|Homo': [50, 335], 'new|TRBV5-7*01_VUEJ|Homo': [50, 335], 'new|TRBV5-8*01_7EBE|Homo': [50, 335],
    'new|TRBV5-8*01_AK3E|Homo': [50, 335], 'new|TRBV5-8*01_H3GH|Homo': [50, 335], 'new|TRBV5-8*01_RS7C|Homo': [50, 335], 'new|TRBV5-8*01_YC2B|Homo': [50, 335],
    'new|TRBV6-2*01_AM56|Homo': [50, 336], 'new|TRBV6-2*01_CBCI|Homo': [50, 336], 'new|TRBV6-2*01_HJ2Z|Homo': [50, 336], 'new|TRBV6-2*01_XTOS|Homo': [50, 335],
    'new|TRBV6-4*01_GP7T|Homo': [50, 336], 'new|TRBV6-5*01_TJKV|Homo': [50, 328], 'new|TRBV6-6*01_UD2U|Homo': [50, 336], 'new|TRBV6-6*01_FIVH|Homo': [50, 336],
    'new|TRBV6-7*01_7LMT|Homo': [50, 336], 'new|TRBV6-8*01_KTY7|Homo': [50, 333], 'new|TRBV6-8*01_XLO6|Homo': [50, 333], 'new|TRBV6-9*01_EZAS|Homo': [50, 336],
    'new|TRBV6-9*01_P574|Homo': [50, 336], 'new|TRBV7-1*01_P2EZ|Homo': [50, 339], 'new|TRBV7-3*01_45HE|Homo': [50, 339], 'new|TRBV7-3*01_CONF|Homo': [50, 339],
    'new|TRBV7-4*01_5GOR|Homo': [50, 339], 'new|TRBV7-4*01_FWAD|Homo': [50, 339], 'new|TRBV7-4*01_QJUC|Homo': [50, 339], 'new|TRBV7-5*02_D26P|Homo': [50, 337],
    'new|TRBV7-5*02_H5DL|Homo': [50, 337], 'new|TRBV7-5*02_KKBI|Homo': [50, 337], 'new|TRBV7-6*01_DNDZ|Homo': [50, 339], 'new|TRBV7-6*01_LZ2G|Homo': [50, 339],
    'new|TRBV7-7*01_PU7O|Homo': [50, 339], 'new|TRBV7-7*01_VURM|Homo': [50, 339], 'new|TRBV7-7*01_4BVP|Homo': [50, 339], 'new|TRBV7-8*01_O2BF|Homo': [50, 339],
    'new|TRBV7-9*03_2FBC|Homo': [50, 339], 'new|TRBV7-9*03_BU7F|Homo': [50, 339], 'new|TRBV7-9*03_HU76|Homo': [50, 339], 'new|TRBV7-9*03_MUA6|Homo': [50, 339],
    'new|TRBV8-1*02_OOYH|Homo': [50, 328], 'new|TRBV8-2*02_VBY5|Homo': [50, 322], 'new|TRBVA*01_7YST|Homo': [50, 311], 'new|TRBVA*02_PCME|Homo': [50, 311],
    'new|TRBVB*01_H7AX|Homo': [50, 401], 'new|TRBVB*02_I3CD|Homo': [50, 401], 'new|TRBVB*02_Z2VV|Homo': [50, 401], 'new|TRBVC*01_O7LI|Homo': [50, 138],
    'new|TRDV1*01_LOTC|Homo': [50, 336], 'new|TRDV3*01_HHZX|Homo': [50, 339], 'new|TRGV10*02_GFCU|Homo': [50, 358], 'new|TRGV10*02_NFML|Homo': [50, 358],
    'new|TRGV10*02_QIWQ|Homo': [50, 358], 'new|TRGV10*02_FO2S|Homo': [50, 359], 'new|TRGV11*01_PB6P|Homo': [50, 358], 'new|TRGV11*01_ZDXH|Homo': [50, 358],
    'new|TRGV11*01_RA2E|Homo': [50, 364], 'new|TRGV2*01_G6VX|Homo': [50, 349], 'new|TRGV3*01_5AGX|Homo': [50, 349], 'new|TRGV3*01_623Z|Homo': [50, 349],
    'new|TRGV3*01_EAGP|Homo': [50, 349], 'new|TRGV3*01_XARG|Homo': [50, 349], 'new|TRGV3*01_ZZNT|Homo': [50, 349], 'new|TRGV4*02_5ZCB|Homo': [50, 349],
    'new|TRGV4*02_BLOG|Homo': [50, 349], 'new|TRGV4*02_C5MV|Homo': [50, 349], 'new|TRGV4*02_CZIJ|Homo': [50, 349], 'new|TRGV4*02_K52A|Homo': [50, 349],
    'new|TRGV5*01_EVVF|Homo': [50, 349], 'new|TRGV5P*02_6TPL|Homo': [50, 349], 'new|TRGV5P*02_E3DS|Homo': [50, 349], 'new|TRGV5P*02_IZMC|Homo': [50, 349],
    'new|TRGV6*01_IJXA|Homo': [50, 348], 'new|TRGV6*01_TUMS|Homo': [50, 348], 'new|TRGV7*01_Y3NC|Homo': [50, 348], 'new|TRGV8*01_3G66|Homo': [50, 349],
    'new|TRGV9*01_7TBT|Homo': [50, 355], 'new|TRGV9*01_SKJQ|Homo': [50, 355], 'new|TRGVA*01_3DSM|Homo': [50, 334], 'new|TRGVA*01_4MG3|Homo': [50, 334],
    'new|TRGVA*01_A6JA|Homo': [50, 334], 'new|TRGVA*01_N2NN|Homo': [50, 334], 'new|TRGVB*01_24EP|Homo': [50, 362], 'new|TRGVB*02_A25Y|Homo': [50, 361],
    'new|TRGVB*02_X42H|Homo': [50, 361], 'new240314|TRAV1-1*01_FYFB|Homo': [50, 324], 'new240314|TRAV1-1*01_YZ7J|Homo': [50, 324], 'new240314|TRAV12-2*01_XEQC|Homo': [50, 326],
    'new240314|TRAV12-2*01_OPDL|Homo': [50, 326], 'new240314|TRAV14/DV4*01_7OPH|Homo': [50, 339], 'new240314|TRAV14/DV4*01_J3MI|Homo': [50, 339], 'new240314|TRAV14/DV4*01_H3S5|Homo': [50, 339],
    'new240314|TRAV21*01_2KLY|Homo': [50, 328], 'new240314|TRAV35*01_QREI|Homo': [50, 325], 'new240314|TRAV35*01_UOVN|Homo': [50, 325], 'new240314|TRAV35*01_4ULT|Homo': [50, 325],
    'new240314|TRAV35*01_BZFC|Homo': [50, 325], 'new240314|TRAV36/DV7*05_HFAO|Homo': [50, 326], 'new240314|TRAV36/DV7*05_CAO3|Homo': [50, 326], 'new240314|TRAV38-1*01_4K2K|Homo': [50, 339],
    'new240314|TRAV38-1*01_WAQH|Homo': [50, 339], 'new240314|TRAV6*01_N2UE|Homo': [50, 329], 'new240314|TRAV6*01_PCAP|Homo': [50, 329], 'new240314|TRAV6*01_ZQU6|Homo': [50, 329],
    'new240314|TRAV8-4*01_2KEZ|Homo': [50, 333], 'new240314|TRAV8-4*01_MEVB|Homo': [50, 333], 'new240314|TRAV9-2*02_K7YF|Homo': [50, 330], 'new240314|TRAV9-2*02_QXIT|Homo': [50, 330],
    'new240314|TRBV11-3*01_JDJJ|Homo': [50, 339], 'new240314|TRBV19*01_V7QQ|Homo': [50, 336], 'new240314|TRBV24/OR9-2*03_6MOS|Homo': [50, 342], 'new240314|TRBV5-5*01_KIIP|Homo': [50, 335],
    'new240314|TRBV5-5*01_EABI|Homo': [50, 335], 'new240314|TRBV7-5*02_3DWO|Homo': [50, 338], 'new240314|TRGV2*01_DZZ3|Homo': [50, 349], 'new240314|TRGV3*01_GTOV|Homo': [50, 349],
    'new240712|TRAV6*01_KBZX|Homo': [50, 329], 'new240712|TRBV10-1*02_HA27|Homo': [50, 336], 'new240712|TRBV20-1*01_SXO6|Homo': [50, 342], 'new240712|TRBV20/OR9-2*03_V4XL|Homo': [50, 342],
    'new240712|TRBV10-3*01_RK3P|Homo': [50, 336], 'new240712|TRBV14*01_OB3I|Homo': [50, 339], 'new240712|TRBV20-1*01_V6LC|Homo': [50, 342], 'new240712|TRBV5-4*01_JP3T|Homo': [50, 335],
    'new240712|TRDV2*01_K5BK|Homo': [50, 337], 'new240712|TRBV20/OR9-2*03_3HNI|Homo': [50, 342], 'new240712|TRAV23/DV6*01_TGPA|Homo': [50, 329], 'new240712|TRBV4-3*01_3IWL|Homo': [50, 336],
    'new240712|TRAV30*01_S4SX|Homo': [50, 323], 'new240712|TRBV20-1*01_WDZ6|Homo': [50, 342], 'new240712|TRBV30*02_6AST|Homo': [50, 333], 'new240712|TRBV20/OR9-2*03_XPQ2|Homo': [50, 342],
    'new240712|TRBV5-4*01_MMDA|Homo': [50, 335], 'new240712|TRBV6-6*01_UXBN|Homo': [50, 336], 'new240712|TRBV30*02_OTVF|Homo': [50, 333], 'new240712|TRBV5-4*01_23O6|Homo': [50, 335],
    'new240712|TRBV30*02_73RW|Homo': [50, 333]}
    
    return FL_mark_dict


from difflib import SequenceMatcher

def assign_new_name_basic(basic_name, dict_target):
    ext_num = 0
    while True:
        tmp_name = basic_name + '_' + str(ext_num)
        if dict_target.get(tmp_name):
            ext_num+=1
        else:
            return tmp_name

def parse_fasta2(fn_fasta):
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
    '''parse the fasta file into a dictionary'''
    # dict_name_SEQ {}
    #  - keys: seq_name
    #  - values: seq_SEQ
    dict_name_SEQ = {}
    dict_novel_only_SEQ = {}

    def add_sequence(seq_name, seq_SEQ):
        if seq_name in dict_name_SEQ:
            print("WARNING! Duplicate sequence name:", seq_name)
            seq_name = assign_new_name_basic(seq_name, dict_name_SEQ)
        dict_name_SEQ[seq_name] = seq_SEQ

        # Check if the sequence is novel
        if "novel" in seq_name:
            dict_novel_only_SEQ[seq_name] = seq_SEQ

    with open(fn_fasta, 'r') as f_f:
        seq_name = ""
        seq_SEQ = ""
        
        for line in f_f:
            line = line.strip()
            if line.startswith('>'):
                if seq_name:
                    add_sequence(seq_name, seq_SEQ)
                seq_name = line[1:].split()[0]
                seq_SEQ = ""
            else:
                seq_SEQ += line
        
        if seq_name: # Add the last sequence if it exists
            add_sequence(seq_name, seq_SEQ)

    return dict_name_SEQ, dict_novel_only_SEQ


def record_differences(sequence1, sequence2):
    differences = []
    for i, (char1, char2) in enumerate(zip(sequence1, sequence2)):
        if char1 != char2:
            differences.append(i)
    return differences



def find_closest_match(query_sequence, fasta_dict, query_gene):
    '''find the closest match for the {input query}'''
    closest_matches = []
    closest_match_ratios = []
    # query_seq vs. subject_seq
    for key, sequence in fasta_dict.items():
        if key != query_gene:  # skip the query
            matcher = SequenceMatcher(autojunk=False, a=query_sequence.lower(), b=sequence.lower())
            match_ratio = matcher.ratio()
            if match_ratio >= 0.5: # if the match ratio is greater than 0.5, add the gene name and match ratio to the list
                closest_matches.append(key)
                closest_match_ratios.append(match_ratio)

    return closest_matches


def compare_sequences(fasta_dict):
    '''# key: gene name; # value: list of positions where the sequences differ'''
    differences = {} #record only the positions where the sequences differ
    differences_full = {} #record positions and names of subjects (closest matches)
    for query_gene, query_sequence in fasta_dict.items(): # key: gene name; sequence: sequence of the gene
        if "novel" in query_gene:
            break
        ## compare the sequence with all other sequences
        ### closest matches: list of closed matched (>=50%) allele names
        closest_matches = find_closest_match(query_sequence, fasta_dict, query_gene) 
        
        ## to record the differences between the query seq and close subject seq list
        ### key: query_gene
        ### value: list of positions where the sequences differ
        differences[query_gene] = []
        differences_full[query_gene]=[]
        
        for close_allele in closest_matches:
            subject_sequence = fasta_dict[close_allele] # sequence2: sequence of every closest match (seq of every allele)
            # querySeq_aligned, subjectSeq_aligned = align_sequences(query_sequence, subject_sequence) # align query sequence and subject sequence for comparison
            # diff_positions = record_differences(querySeq_aligned, subjectSeq_aligned) # diff_positions: list of positions where the sequences differ
            diff_positions = record_differences(query_sequence, subject_sequence) # diff_positions: list of positions where the sequences differ
            diff_positions = [pos for pos in diff_positions if (pos <= (len(query_sequence)-1)) and (len(diff_positions)<10)] # remove positions that are beyond the original length (in order to record the differences on the query sequence only)

            ### differences[query_gene]: [1,2,3...] (list of positions where the sequences differ)            
            differences[query_gene].append(diff_positions) # record the different position
            differences_full[query_gene].append(f'vs. {close_allele}: {diff_positions}') # record the different position and the name of the closest match           
    for query_gene in differences:
        differences[query_gene] = sorted(list(set([element for sublist in differences[query_gene] for element in sublist])))
        print(f'>>> {query_gene} : {differences_full[query_gene]}')
    
    return differences


def compare_genes(fasta_dict, new_fasta_dict):
    '''merge known and novel fasta dictionaries'''
    for new_gene in new_fasta_dict:
        if new_gene not in fasta_dict:
            fasta_dict[new_gene] = new_fasta_dict[new_gene]
    
    '''compare sequences and return differences
    # key: gene name
    # value: list of positions where the sequences differ
    '''
    allelePos_dict= compare_sequences(fasta_dict)


    return allelePos_dict
