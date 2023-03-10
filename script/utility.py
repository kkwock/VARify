# libraries
import warnings
import pandas as pd


# For every row in here, we're going to apply this uniqueSNP function
def mergeTab(mp, so):
    sub_mp = mp[['snp_pos', 'pileup']]
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


def uniqueSNP(x):
    l = str(x['pileup'])
    # make upper
    l = l.upper()  # just to make simple
    l = varOnly(l)

    # Get list of unique
    a = [i for i in l]
    a_unique = list(set(a))

    # create a dictionary
    df = pd.DataFrame()
    temp = [i for i in x['ref_codon']]
    snp_dict = {}

    flag = False

    if len(a_unique) > 2:
        warnings.warn("There are more than 2 SNPS at position");
        print(x['snp_pos'])
        flag = True

    for snp in a_unique:
        snp_dict['snp_pos'] = int(x['snp_pos'])
        snp_dict['varify_alt'] = snp

        if flag:
            snp_dict['flag'] = '*'
        else:
            snp_dict['flag'] = ''

        for i in range(len(temp)):
            temp[i] = snp
            snp_dict[f"varify_codon{i + 1}"] = ''.join(temp)
            temp = [i for i in x['ref_codon']]

        df = df.append(snp_dict, ignore_index=True)

        snp_dict = {}

    df = pd.merge(sub_cf, df, on='snp_pos')
    return df


def varify(x):
    df = x.apply(axis=1, func=uniqueSNP)
    df = pd.concat(list(df))
    reorg = ['chr_id', 'snp_pos', 'ref_allele', 'alt_allele',
             'varify_alt', 'codon1_genome_pos', 'codon2_genome_pos',
             'codon3_genome_pos', 'ref_codon', 'varify_codon1',
             'varify_codon2', 'varify_codon3', 'flag']

    return df[reorg]