import argparse as ap
import pysam
import sys

def parse_args():
    parser = ap.ArgumentParser(prog="MIRCHI")
    parser._action_groups.pop()
    print("MIRCHI - Mendelian Inheritance Repeat Concordance in Human Individuals\n")
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-k', '--kid', required=True, metavar='<FILE>', help='offspring VCF file.')
    required.add_argument('-d', '--dad', required=True, metavar='<FILE>', help='paternal VCF file.')
    required.add_argument('-m', '--mom', required=True, metavar='<FILE>', help='maternal VCF file.')
    required.add_argument('-r', '--regions', required=True, metavar='<FILE>', help='regions BED file in ATaRVa format.')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--non-ref', action='store_true', help='If set, will not consider "same as reference alleles" in the analysis. [default: All-loci-mode]')
    optional.add_argument('-o', '--out', metavar='<FILE>', default='mirchi.tsv', help='output file name.')
    optional.add_argument('--contigs', nargs='+', help='contigs to calculate mendelian concordance eg [chr1 chr2 ch3]. If not mentioned every contigs in the region file will be considered.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()


def extract_motif(tbx, contigs, non_ref):
    Assembly_motif_dict = {}
    atarva_non_ref = {}
    for each_chr in contigs:
        for line in tbx.fetch(each_chr):
            line = line.strip().split('\t')
            key = f'{line[0]}:{line[1]}-{line[2]}'
            atarva_non_ref[key] = [not non_ref, False, False]
            Assembly_motif_dict[key] = int(float(line[4]))
    return Assembly_motif_dict, atarva_non_ref

def extract_alleles(vcf, contig, atarva_non_ref, allele_dict):

    for line in vcf.fetch(contig):
        line = line.strip().split('\t')
        if line[6]!='PASS': continue
        start = int(line[1])
        end = int(line[7].split(';')[3].split('=')[1])
        key = f'{line[0]}:{start}-{end}'
        sample = line[-1].split(':')
        alleles = sample[1].split(',')
        allele_set = set()
        for i in alleles:
            allele_set.add(int(i))
        if (len(allele_set)==1) and (end-start == list(allele_set)[0]):
            pass
        else:
            atarva_non_ref[key][0] = True
        allele_dict[key] = allele_set
    return allele_dict  

def checker(motif, diff_in_F, diff_in_M):
    men_al = 0
    one_al = 0
    mot_al = 0
    no_men = 0
    
    if diff_in_F == 0:
        men_al += 1
    elif diff_in_F == 1:
        one_al += 1
    elif diff_in_F <= motif:
        mot_al += 1
    else:
        no_men+=1

    if diff_in_M == 0:
        men_al += 1
    elif diff_in_M == 1:
        one_al += 1
    elif diff_in_M <= motif:
        mot_al += 1
    else:
        no_men+=1

    return [men_al, one_al, mot_al, no_men]

def mendelian_concordance(A_Whole_loci_set, atarva_non_ref, dad_alleles, mom_alleles, kid_alleles, Assembly_motif_dict):
    total_alleles = 0
    mendelian_alleles = 0
    non_mendelian_alleles = 0
    common_loci = 0

    off_by_one = 0
    off_by_1motif = 0 

    for each_loci in A_Whole_loci_set:
        if any(atarva_non_ref[each_loci]) and (each_loci in dad_alleles) and (each_loci in mom_alleles) and (each_loci in kid_alleles):

            motif = Assembly_motif_dict[each_loci] # gets motif size from already generated dictionary
            common_loci += 1
            total_alleles += 2 # for each locus, total alleles is incremented by 2

            son_alleles = list(kid_alleles[each_loci])
            if len(son_alleles)==1: son_alleles=son_alleles*2
            father_alleles = list(dad_alleles[each_loci])
            if len(father_alleles)==1: father_alleles=father_alleles*2
            mother_alleles = list(mom_alleles[each_loci])
            if len(mother_alleles)==1: mother_alleles=mother_alleles*2

            O1 = son_alleles[0] # son allele 1
            O2 = son_alleles[1] # son allele 2
            # calculating the difference between son allele and parent allele in two combination
            case1_diff = ([abs(O1-i) for i in father_alleles], [abs(O2-i) for i in mother_alleles])
            case2_diff = ( [abs(O2-i) for i in father_alleles], [abs(O1-i) for i in mother_alleles])
            # Summation of the min values from each combination
            O1f_O2m = sum([min([abs(O1-i) for i in father_alleles]), min([abs(O2-i) for i in mother_alleles])])
            O1m_O2f = sum([min([abs(O2-i) for i in father_alleles]), min([abs(O1-i) for i in mother_alleles])])
            
            sum_diff_pair = [O1f_O2m, O1m_O2f] # list of (sum of min diff values)
            if O1f_O2m!=O1m_O2f:  # if the values are not same, take the combination which gives the minimum value
                idx = sum_diff_pair.index(min(sum_diff_pair)) # getting the index of min value of sum diff
                min_pair = sum_diff_pair[idx] # getting the min value
                if min_pair == 0: # if the min value is 'Zero', its a perfect mendelian alleles
                    mendelian_alleles += 2
                elif idx==0: # if the min value is from 1st combo, check their class and increment it
                    diff_in_F = min(case1_diff[0])
                    diff_in_M = min(case1_diff[1])
                    cat_al = checker(motif, diff_in_F, diff_in_M)
                    mendelian_alleles += cat_al[0]
                    off_by_one += cat_al[1]
                    off_by_1motif += cat_al[2]
                    non_mendelian_alleles += cat_al[3]
                elif idx==1: # if the min value is from 2nd combo, check their class and increment it
                    diff_in_F = min(case2_diff[0])
                    diff_in_M = min(case2_diff[1])
                    cat_al = checker(motif, diff_in_F, diff_in_M)
                    mendelian_alleles += cat_al[0]
                    off_by_one += cat_al[1]
                    off_by_1motif += cat_al[2]
                    non_mendelian_alleles += cat_al[3]
            else: # if the two (sum of min diff)s are same
                
                # eg : son_alleles = [14,15]
                #      father_alleles = [14, 17]
                #      mother_alleles = [14, 16]
                
                # check which combo has minimum difference with any of the son allele
                c1 = min([min([abs(O1-i) for i in father_alleles]), min([abs(O2-i) for i in mother_alleles])]) # eg: [0,3]  [1,1], # this line of code gives 0       
                c2 = min([min([abs(O2-i) for i in father_alleles]), min([abs(O1-i) for i in mother_alleles])]) # eg: [0,2]  [1,2], # this line of code gives 0     
                
                if c1<c2: # if the value from 1st combo is smaller, that combo is taken for mendelian classification # eg c1 = 0, c2 = 1
                    diff_in_F = min(case1_diff[0])
                    diff_in_M = min(case1_diff[1])
                    cat_al = checker(motif, diff_in_F, diff_in_M)
                    mendelian_alleles += cat_al[0]
                    off_by_one += cat_al[1]
                    off_by_1motif += cat_al[2]
                    non_mendelian_alleles += cat_al[3]
                elif c2<c1: # if the value from 2nd combo is smaller, that combo is taken for mendelian classification # eg c1 = 1, c1 = 0
                    diff_in_F = min(case2_diff[0])
                    diff_in_M = min(case2_diff[1])
                    cat_al = checker(motif, diff_in_F, diff_in_M)
                    mendelian_alleles += cat_al[0]
                    off_by_one += cat_al[1]
                    off_by_1motif += cat_al[2]
                    non_mendelian_alleles += cat_al[3]
                elif c1==c2: # if both values are same then,
                    # check the 2nd minimum value from both combo
                    m1 = max([min([abs(O1-i) for i in father_alleles]), min([abs(O2-i) for i in mother_alleles])]) # eg: [0,3]  [1,1], # this line of code gives 1 
                    m2 = max([min([abs(O2-i) for i in father_alleles]), min([abs(O1-i) for i in mother_alleles])]) # eg: [0,2]  [1,2], # this line of code gives 1 
                    if m1<m2: # if the value from 1st combo is smaller, that combo is taken for mendelian classification # eg m1 = 0, m2 = 1
                        diff_in_F = min(case1_diff[0])
                        diff_in_M = min(case1_diff[1])
                        cat_al = checker(motif, diff_in_F, diff_in_M)
                        mendelian_alleles += cat_al[0]
                        off_by_one += cat_al[1]
                        off_by_1motif += cat_al[2]
                        non_mendelian_alleles += cat_al[3]
                    elif m2<m1: # if the value from 1st combo is smaller, that combo is taken for mendelian classification # eg m1 = 1, m2 = 0
                        diff_in_F = min(case2_diff[0])
                        diff_in_M = min(case2_diff[1])
                        cat_al = checker(motif, diff_in_F, diff_in_M)
                        mendelian_alleles += cat_al[0]
                        off_by_one += cat_al[1]
                        off_by_1motif += cat_al[2]
                        non_mendelian_alleles += cat_al[3]
                    else: # If both values are same, take any one of the combo
                        diff_in_F = min(case1_diff[0])
                        diff_in_M = min(case2_diff[1])
                        cat_al = checker(motif, diff_in_F, diff_in_M)
                        mendelian_alleles += cat_al[0]
                        off_by_one += cat_al[1]
                        off_by_1motif += cat_al[2] 
                        non_mendelian_alleles += cat_al[3]

    return total_alleles, mendelian_alleles, non_mendelian_alleles, common_loci, off_by_one, off_by_1motif

def main():
    args = parse_args()

    # Load the VCF files
    kid_vcf = pysam.TabixFile(args.kid)
    dad_vcf = pysam.TabixFile(args.dad)
    mom_vcf = pysam.TabixFile(args.mom)

    # Load the regions BED file
    tbx  = pysam.Tabixfile(args.regions)
    if not args.contigs:
        contigs = sorted(tbx.contigs)
    else:
        contigs = sorted(args.contigs)
    Assembly_motif_dict, atarva_non_ref = extract_motif(tbx, contigs, args.non_ref)
    tbx.close()

    # Extract alleles from the VCF files
    kid_alleles = {}
    dad_alleles = {}
    mom_alleles = {}
    for contig in contigs:
        dad_alleles = extract_alleles(dad_vcf, contig, atarva_non_ref, dad_alleles)
        mom_alleles = extract_alleles(mom_vcf, contig, atarva_non_ref, mom_alleles)
        kid_alleles = extract_alleles(kid_vcf, contig, atarva_non_ref, kid_alleles)
    kid_vcf.close()
    dad_vcf.close()
    mom_vcf.close()

    A_Whole_loci_set = (dad_alleles.keys() & mom_alleles.keys()) & kid_alleles.keys()
    # Process the VCF files and regions
    total_alleles, mendelian_alleles, non_mendelian_alleles, common_loci, off_by_one, off_by_1motif = mendelian_concordance(A_Whole_loci_set, atarva_non_ref, dad_alleles, mom_alleles, kid_alleles, Assembly_motif_dict)

    # Output results
    with open(args.out, 'w') as out_file:
        out_file.write('#MIRCHI - Mendelian Inheritance Repeat Concordance in Human Individuals\n')
        out_file.write(f"##command=MIRCHI {' '.join(sys.argv)}\n")
        out_file.write('\t'.join(['#Common loci', 'Total alleles', 'Mendelian alleles', 'off by one', 'off by 1motif', 'Non mendelian alleles']) + '\n')
        out_file.write('\t'.join(map(str, [common_loci, total_alleles, mendelian_alleles, off_by_one, off_by_1motif, non_mendelian_alleles])) + '\n')

if __name__ == '__main__':
    main()