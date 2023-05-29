# Utility File

# libraries
import warnings
import pandas as pd

asciiDict = {chr(i): (i - 33) for i in range(33, 74)}
ntRevDict = {"A": "T",
             "T": "A",
             "G": "C",
             "C": "G"}


def snpCheck(strand, snp):
    '''
    Checks for pos/neg strand
    returns SNP
    '''

    if strand == '-':
        return ntRevDict.get(snp)
    else:
        return snp


def varify(x, mp):
    # dictionary of alternate codons from pileup
    # alt_dict = altDict(mp)

    # running the script
    snp_df = x.apply(axis=1, func=snpDict)
    snp_df = pd.concat(list(snp_df))
    df = pd.merge(x, snp_df, on='snp_pos')

    # Add the dictionary
    df = df.apply(axis=1, func=translate)
    df = df.apply(axis=1, func=altValidate)
    df = df.apply(axis=1, func=codonVarify)
    df = df.apply(axis=1, func=aaVarify)
    df = df.drop(columns='pileup', axis=1)

    reorg = ['chr_id', 'snp_pos', 'ref_allele', 'alt_allele', 'gene_id', 'mrna_id',
             'prot_id', 'strand', 'effect', 'snp_cds_pos', 'codon1_genome_pos',
             'codon2_genome_pos', 'codon3_genome_pos', 'snp_aa_pos', 'ref_codon',
             'alt_codon', 'ref_aa', 'alt_aa',  # 'ref',
             'pos1_pileup', 'pos2_pileup', 'pos3_pileup',
             'varify_codon', 'varify_allele',
             'varify_aa', 'nt_VARified', 'codon_VARified', 'aa_VARified', 'comment']
    return df[reorg]


def snpDict(x):
    # import itertools # for permutations
    # print(x.snp_pos)
    # Takes the current position of SNP
    p3 = int(x.codon3_genome_pos)
    p2 = int(x.codon2_genome_pos)
    p1 = int(x.codon1_genome_pos)
    subDict = {p1: '.', p2: '.', p3: '.'}

    # Fills dictionary with SNP of position from mpileup
    res_dict = {key: alt_dict.get(key, subDict[key]) for key in subDict}
    res = ''.join(res_dict.values())

    pos_check = False
    pos_list = [p1, p2, p3]
    #     print(pos_list)
    #     print(subDict)

    # If there are any differences between RNA & DNA, return position
    if pos != False:
        for p in pos_list:
            if p in pos:
                # only add p into pos_dict
                pos_check = True
                # snp = snpCheck(x.strand, alt_dict.get(p))
                pos_dict = {p: alt_dict.get(p)}
                # print(pos_dict)
                # pos_dict = {key:val for key, val in pos_dict.items() if val != '.'}

    # starts from reference
    ref = list(x['ref_codon'])

    # create a new dictionary
    snp_dict = {}
    multi_check = False
    comment = ''
    flag = ''

    # create pileup in df
    for i in range(len(pos_list)):
        cname = ['pos1_pileup', 'pos2_pileup', 'pos3_pileup']

        try:
            snp_dict[cname[i]] = mp[mp.snp_pos == pos_list[i]].pileup.tolist()[0]
        except:
            # if the pileup does not have info, replace with *
            snp_dict[cname[i]] = "*"

    # MultiFlag for complex codons
    if (list(res)).count(".") < 2:
        multi_check = True
        if pos_check == False:
            # if < 2, then there are multi snps
            flag = 'complex_codon'
        elif pos_check == True:
            flag = f'NA'
            comment = f'SNP(s) {pos_dict} only in RNA'

    if multi_check == True:
        # if multi_check true, then complex_codon
        snp_dict['varify_codon'] = flag
        snp_dict['snp_pos'] = x['snp_pos']
        current_pos = int(x['snp_pos'])
        snp_dict['varify_allele'] = alt_dict[int(current_pos)]  # fill with current pos snp
        snp_dict['comment'] = comment

    else:
        current_pos = int(x['snp_pos'])
        idx = list(subDict.keys()).index(current_pos)
        # get ref
        try:
            ref[idx] = snpCheck(x.strand, alt_dict[current_pos])
            snp_dict['varify_allele'] = alt_dict[current_pos]  # fill with current pos snp
        except:
            ref[idx] = snpCheck(x.strand, x.ref_allele)
            snp_dict['varify_allele'] = x.ref_allele

        snp_dict['varify_codon'] = ''.join(ref)
        snp_dict['snp_pos'] = x['snp_pos']
        snp_dict['comment'] = comment

    snp_df = pd.DataFrame([snp_dict])

    return snp_df


def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)


def getAlt(x):
    indel = has_numbers(x.pileup)
    l = list(x['pileup'])
    nts = ['A', 'C', 'T', 'G']

    # check if indel, then skip
    while indel == True:
        for i in l:
            if i.isnumeric():
                string = ''.join(l)
                n = int(i)
                idx = l.index(i)
                r = range(idx - 1, idx + n + 1)
                pat = ''.join(l[r[0]:r[-1] + 1])
                string = string.replace(pat, '')
        l = [i for i in string]
        indel = has_numbers(string)

    l = [i for i in l if i.isalpha()]
    l = [i for i in l if i.upper() in nts]

    if len(l) == 0:
        x['varify_allele'] = None
    else:
        # Get list of unique
        a = [i for i in l]
        a = pd.DataFrame(a)
        # a_unique = list(set(a))

        # getting max alt
        # a_unique = pd.DataFrame(a_unique)

        # Get the max snp
        # nt = a_unique.groupby([0]).apply(lambda x: x.value_counts().index[0])[0]
        # nt = nt[0]
        # nt = pd.value_counts(a_unique[0].values.flatten()).index[0]
        nt = pd.value_counts(a.values.flatten()).index[0]
        x['varify_allele'] = nt.upper()

    return x


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def checkRNADict(x, mp):
    x_dict = altDict(x)  # for DNA
    mp_dict = altDict(mp)  # for RNA

    if str(set(mp_dict) - set(x_dict)) == 'set()':
        return False, mp_dict
    else:
        # returns position difference from RNA
        return list(set(mp_dict) - set(x_dict)), mp_dict


def readMP(x):
    import pandas as pd
    mp = pd.read_table(x, header=None)
    mp = mp.rename(columns={0: "chr", 1: "snp_pos", 2: "ref", 3: "reads", 4: "pileup"})
    mp = mp[mp.reads >= 1]
    mp = mp.dropna(axis=0)
    return (mp)


def aaVarify(x):
    x['aa_VARified'] = 'Yes' if x['varify_aa'] == x['alt_aa'] else "No"
    return x


def codonVarify(x):
    x['codon_VARified'] = 'Yes' if x['varify_codon'] == x['alt_codon'] else 'No'
    return x


def translate(x):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W', }

    codon = x['varify_codon']

    if codon not in table.keys():
        x['varify_aa'] = ""
    else:
        x['varify_aa'] = table[codon]

    return (x)


# For every row in here, we're going to apply this uniqueSNP function

def mergeTab(mp, so):
    sub_mp = mp[['snp_pos', 'pileup']]
    # sub_mp['snp_pos'] = str(sub_mp['snp_pos'])
    mdf = pd.merge(so, sub_mp, on='snp_pos')
    return mdf


def varOnly(a):
    '''
    purpose: remove non-nt values in pileup
    input: string of mpileup
    output: nt string
    '''
    x = [i for i in a if i.upper() in ['A', 'C', 'T', 'G']]

    return x


def isMono(ref, temp):
    # Check if two snps are monomorphic
    c = 0
    count = [c + 1 for snp in temp if snp == ref]
    if sum(count) >= 2:
        warnings.warn("More than 2 monomorphic alleles at position.")
        return True
    else:
        return False


def altValidate(x):
    # Varify alt == Alt? T or F
    # Add flag to capture

    x['nt_VARified'] = 'Yes' if x['alt_allele'] == x['varify_allele'] else 'No'

    return x


def effect(x):
    # missense
    # nonsense
    x['effect'] = '_'


###################Older, probably outdated############################

# Codon Combos
import numpy as np
import math


def altDict(x):
    '''
    Takes mPileup
    Out: mpileup dictionary with snp_pos and varify_allele
    '''

    # dataframe of mpileup goes through each pileup with apply
    x = x.apply(axis=1, func=getAlt)
    alt_dict = dict(zip((x['snp_pos']), x['varify_allele']))
    clean_dict = {k: alt_dict[k] for k in alt_dict if pd.isna(alt_dict[k]) == False}

    return clean_dict


def flatten(d, parent_key='', sep='_'):
    import collections
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

# import itertools

# permutations = list(itertools.product(["A","G"], ["T","C"], ["A"]))

# print(permutations)

# ["".join(x) for x in permutations]