import pandas as pd
import re
import argparse
import json
from datetime import datetime
import concurrent.futures
import multiprocessing
import pickle


def sam2pair_extract(s2p_out):
    # Extract every 1st 2nd and 4th entry into separate DataFrames
    col1 = s2p_out.iloc[::4].reset_index(drop=True)
    col2 = s2p_out.iloc[1::4].reset_index(drop=True)
    col4 = s2p_out.iloc[3::4].reset_index(drop=True)

    # Split the data in column 1 by tab delimiter and keep elements at indices 1, 3, 4, and 9
    split_col1 = col1[0].str.split('\t', expand=True).iloc[:, [0, 1, 2, 3, 4, 5]]

    # Getting softclip number
    split_col1[5] = split_col1[5].apply(softclip_number)

    # Create a new DataFrame with the split data and col2_cleaned as columns
    result_df = pd.concat([split_col1, col2, col4], axis=1)

    # Rename the columns
    result_df.columns = ['QNAME', 'FLAG', 'CHR', 'POS', 'MAPQ', 'SOFTCLIP', 'SEQ', 'REF']

    result_df = result_df.dropna(subset=['QNAME'])

    # remove 'S' from REF ()
    result_df['REF'] = result_df.apply(lambda row: remove_softclip(row, 'REF'), axis=1)
    result_df['SEQ'] = result_df.apply(lambda row: remove_softclip(row, 'SEQ'), axis=1)

    # Get Range
    result_df['POS_END'] = result_df.apply(lambda row: int(row['POS']) + len(row['SEQ']), axis=1)

    # Interpret the SAM flag
    result_df['FLAG_EXPLAINED'] = result_df['FLAG'].astype(int).apply(sam_flag_explainer)

    #Drop unmapped
    result_df = result_df[result_df['FLAG_EXPLAINED'].apply(lambda x: x.get('read_unmapped', False) != True)]


    return result_df[['QNAME', 'FLAG', 'CHR', 'POS', 'POS_END', 'MAPQ', 'SOFTCLIP', 'SEQ', 'REF', 'FLAG_EXPLAINED']]


def sam_flag_explainer(flag:int):
    """
    Intepret sam flag and make a flag-explained-dictionary.

    """
    try:
        # Convert sam flags to binary format, like 0101001
        # flag is a integer.
        flag_binary = "{0:b}".format(flag)
    except Exception as err:
        if int(flag) == flag:
            flag = int(flag)
            flag_binary = "{0:b}".format(flag)
        else:
            print(f"{flag=}")
            raise ValueError from err

    flag_template = {
        "read_paired":False,
        "read_mapped_in_proper_pair":False,
        "read_unmapped":False,
        "mate_unmapped":False,
        "read_reverse_strand":False,
        "mate_reverse_strand":False,
        "first_in_pair":False,
        "second_in_pair":False,
        "not_primary_alignment":False,
        "read_fails_platform/vendor_quality_checks":False,
        "read_is_PCR_or_optical_duplicate":False,
        "supplementary_alignment":False,
    }

    for k, f in zip(flag_template.keys(), flag_binary[::-1]): # should intepret the flag from back.
        if int(f) == 1:
            flag_template[k] = True

    return flag_template


def softclip_number(input_string):
    pattern = r'^(\d+)S'

    # Use re.match to find the pattern at the beginning of the string
    match = re.match(pattern, input_string)

    if match:
        # Extract the matched part (the number before 'S')
        result = match.group(1)
        return result

    else:
        return "0"

def remove_softclip(row, type = ['SEQ', 'REF']):
    s_value = int(row['SOFTCLIP'])
    ref_value = row[type]

    if s_value:
        return ref_value[s_value:]
    else:
        return ref_value

def codon_perc(codons, sort_by_percentage=True):
    perc = {}

    for chr_id, targets in codons.items():
        perc[chr_id] = {}
        for target, codon_info in targets.items():
            perc[chr_id][target] = {'percentage': {}}
            codons_list = list(codon_info['codons'])
            codon_counts = {}

            for c1 in codons_list:
                codon_counts[c1] = codon_counts.get(c1, 0) + 1

            for c2, count in codon_counts.items():
                percentage = (count / sum(codon_counts.values())) * 100
                perc[chr_id][target]['percentage'][c2] = percentage

            if sort_by_percentage:
                perc[chr_id][target]['percentage'] = dict(sorted(perc[chr_id][target]['percentage'].items(), key=lambda x: x[1], reverse=True))

    return perc

def get_codons(row, snpeff_json, codons):
    '''
    1. Go through each s2p record
    2. Check if CHR in chr_id
    3. Check if snp_pos between POS:POS_END
        a. Yes: calculate positions and grab SNP
        b. No: Skip
    '''
    start = int(row['POS']) - 1
    chr_id = row['CHR']

    # Check if chr_id exists in the dictionary
    if chr_id not in codons:
        codons[chr_id] = {}

    # Consolidate the dictionaries into one
    for target in snpeff_json[chr_id]:
        if int(target) in range(int(row['POS']), int(row['POS_END'])):
            pos = int(target) - start

            pos1 = snpeff_json[chr_id][target]["codon1_genome_pos"] - start
            pos2 = snpeff_json[chr_id][target]["codon2_genome_pos"] - start
            pos3 = snpeff_json[chr_id][target]["codon3_genome_pos"] - start

            strand = snpeff_json[chr_id][target]['strand']

            try:
                alt1 = row['SEQ'][get_pos(row, pos1)]
            except IndexError:
                alt1 = "X"

            try:
                alt2 = row['SEQ'][get_pos(row, pos2)]
            except IndexError:
                alt2 = "X"

            try:
                alt3 = row['SEQ'][get_pos(row, pos3)]
            except IndexError:
                alt3 = "X"

            # Check if target exists in the dictionary
            if target not in codons[chr_id]:
                codons[chr_id][target] = {'codons': []}

            # check if positive or negative strand
            ntRevDict = {"A": "T",
             "T": "A",
             "G": "C",
             "C": "G",
             "X": "X",
             "-": "-"}

            if strand == '-':
                codons[chr_id][target]['codons'].append(f"{ntRevDict.get(alt1)}{ntRevDict.get(alt2)}{ntRevDict.get(alt3)}")
            else:
                # Add the codon to the target
                codons[chr_id][target]['codons'].append(f"{alt1}{alt2}{alt3}")


    for chr_id, targets in codons.items():
        for target, codon_info in targets.items():

            codons[chr_id][target]['codons'] = list((codon_info['codons']))


    return codons

def get_pos(x, target):
    index = 0
    letter_count = 0
    for char in x['REF']:
        if char != "-":
            letter_count += 1
        if letter_count == target:
            break
        index += 1

    return index

def varify_codons(snpeff_table, codons):
    chr_id = snpeff_table['chr_id']
    snp_pos = str(snpeff_table['snp_pos'])

    snpeff_table['varify_codons'] = "NA"

    if chr_id in codons.keys():
        if snp_pos in codons[chr_id].keys():
            snpeff_table['varify_codons'] = codons[chr_id][snp_pos]['percentage']
    return snpeff_table


def to_json(snpeff_df):
    json_data = {}

    # Iterate over DataFrame rows and organize data into structured format
    for index, row in snpeff_df.iterrows():
        chr_id = row["chr_id"]
        snp_pos = row["snp_pos"]
        codon1_genome_pos = row["codon1_genome_pos"]
        codon2_genome_pos = row["codon2_genome_pos"]
        codon3_genome_pos = row["codon3_genome_pos"]
        strand = row['strand']

        # Check if the chromosome ID already exists in the structured data
        if chr_id not in json_data:
            json_data[chr_id] = {}

        # Check if the SNP position already exists under the chromosome ID
        if snp_pos not in json_data[chr_id]:
            json_data[chr_id][snp_pos] = {
                "codon1_genome_pos": codon1_genome_pos,
                "codon2_genome_pos": codon2_genome_pos,
                "codon3_genome_pos": codon3_genome_pos,
                'strand':strand
            }

    return json_data

def is_varify(x):
    x['alt_expressed?'] = x.alt_codon in x.varify_codons
    return x

def var_split(x):
    import ast
    # x is the dataframe
    dict_string = x.varify_codons
    alt = x.alt_codon
    dictionary = ast.literal_eval(dict_string)

    # VARify Alt
    varify_alt = {}
    other_alt = {}
    for key, value in dictionary.items():
        if key == alt:
            varify_alt[key] = value
        elif key != alt:
            other_alt[key] = value

    if varify_alt: 
        x['alt_codon_VARified'] = varify_alt
    else: 
        x['alt_codon_VARified'] = "NA"

    if other_alt:
        x['alt_codon_not_VARified'] = other_alt
    else:
        x['alt_codon_not_VARified'] = "NA"

    return x

def main(args):
    date = datetime.now()
    s2pout = args.sam2pairwise
    snpeff = args.snpeff

    # read in data
    s2p = pd.read_csv(s2pout, header=None) # need to ensure None header to work
    snpeff_table = pd.read_csv(snpeff, sep='\t')

    # extract sam2pairwise data
    result_df = sam2pair_extract(s2p)
    result_df.to_csv(f"{date.strftime('%y%m%d_%f')}_sam2pair_table.csv", index=False)

    # Subset snpeff-table data
    sub = ['chr_id', 'snp_pos', 'codon1_genome_pos', 'codon2_genome_pos', 'codon3_genome_pos', 'strand']
    snpeff_df = snpeff_table[sub]
    snpeff_df = snpeff_df.drop_duplicates(subset= sub)

    # Convert snpeff data into JSON for position parsing
    json_data = to_json(snpeff_df)
    json_dump = json.dumps(json_data, indent=4)
    snpeff_json = json.loads(json_dump)

    # Create Codons dictionary
    codons = {}

    for index, row in result_df.iterrows():
        codons = get_codons(row, snpeff_json, codons)

    # add codon percentages to snpeff_table
    perc = codon_perc(codons)
    
    snpeff_table = snpeff_table.apply(varify_codons, codons=perc, axis=1)

    snpeff_table = snpeff_table.apply(is_varify, axis=1)

    # output report
    snpeff_table.to_csv(f"{date.strftime('%y%m%d_%f')}_varify_report.csv", index=False)


if __name__ == '__main__':
    # inputs
    parser = argparse.ArgumentParser(description='VARify Extract Codons')
    parser.add_argument('sam2pairwise', type=str, help='sam2pairwise output')
    parser.add_argument('snpeff', type=str, help='SNPEff-like table')

    # call variables
    args = parser.parse_args()

    main(args)