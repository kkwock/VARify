# Utility File

# libraries
import warnings
import pandas as pd



def aaVarify(x):
    x['aa_VARified'] = 'Yes' if x['varify_aa'] == x['alt_aa'] else "No"
    return x

def codonVarify(x):
    x['codon_VARified'] = 'Yes' if x['varify_codon'] == x['alt_codon'] else 'No'
    return x

def translate(x):
    table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',}
    
    codon = x['varify_codon']

    if codon not in table.keys():
        x['varify_aa'] = ""
    else:
        x['varify_aa'] =table[codon]
        
    return(x)


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
    x = [i for i in a if i.upper() in ['A','C','T','G']]
    
    return x

def isMono(ref, temp):
    # Check if two snps are monomorphic
    c=0
    count = [c+1 for snp in temp if snp == ref]
    if sum(count) >= 2: 
        warnings.warn("More than 2 monomorphic alleles at position.")
        return True
    else:
        return False
    
def multiError(x, pos):
    if len(x) > 2: 
        warnings.warn("There are more than 2 SNPS at position"); print(pos)
    
def uniqueSNP(x):
    '''
    Takes mPileup
    Confirm Alt
    '''
    l = str(x['pileup'])
    # make upper
    l = l.upper() # just to make simple
    l = varOnly(l)

    # Get list of unique
    a = [i for i in l]
    a_unique = list(set(a)) 
    
    # create a dictionary
    df = pd.DataFrame()
    temp = [i for i in x['ref_codon']]
    snp_dict = {}
    
#     flag = False
    
#     multiError(a_unique, x['snp_pos'])

    test = []
    '''
    Creates a combination of SNPs
    Replaces Alt at each position
    '''
    for snp in a_unique:
        snp_dict['snp_pos'] = int(x['snp_pos'])
        snp_dict['varify_allele'] = snp

#         if flag:
#             snp_dict['flag'] = '*'
#         else:
#             snp_dict['flag'] = ''
        
#         for i in range(len(temp)):
#             temp[i] = snp
#             test.append(''.join(temp))
#             snp_dict[f"varify_codon{i+1}"] = ''.join(temp)
#             temp = [i for i in ref]
            
#         flag = isMono(ref, test)  
        
#         if flag: 
#             snp_dict['flag'] = '*'
#         else:
#             snp_dict['flag'] = ''
            
        df = df.append(snp_dict, ignore_index=True)


        snp_dict = {}
        test = []

    df = pd.merge(cf, df, on='snp_pos')
    return df

def altValidate(x):
    # Varify alt == Alt? T or F
    # Add flag to capture
    
    x['nt_VARified'] = 'Yes' if x['alt_allele'] == x['varify_allele'] else 'No'
    
    return x

def varify(x):
    # dictionary of alternate codons from pileup
    alt_dict = altDict(x)

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
    
    return df

def getAlt(x):
    
    '''
    Takes mPileup
    Out: list of uniqueSNP
    '''
    l = x['pileup'].astype(str)[0]
    # make upper
    l = [i for i in l if i.upper() in ['A','C','T','G']]

    # Get list of unique
    a = [i for i in l]
    a_unique = list(set(a)) 
    a_unique
    
    return(a_unique)



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
    
    Creates a combination of SNPs
    if snps have existing alts
    
    Outputs into a dictionary
    
    '''
    # At each line, takes reference codon
    # sub_x = x.loc[x['VARified'] == 'Yes']
    
    # create a dictionary of Alts
    alt_dict = dict(zip((x['snp_pos']), x['alt_allele']))
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