# v.1.0.0

import argparse
import sys
import pandas as pd
import warnings
from pathlib import Path

script_dir = Path(__file__).parent
sys.path.append(script_dir)

import utility as utils

# Input Data
'''
-- in1: mpileup.txt #mpileup table
-- in2: snp_out.txt #SNPEff-like table
'''

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('--in_mp', type=str, required=True)
    PARSER.add_argument('--in_so', type=str, required=True)
    ARGS = PARSER.parse_args()

    # INPUT_CSV = os.path.abspath(ARGS.input_csv)
    # file_ending = INPUT_CSV.split(".")[-1]
    # if file_ending != "txt":
    #     raise Exception("Invalid input file format: File is not a CSV file")
    #
    # csv_input_file = utils.read_file_with_encoding(INPUT_CSV)
    # header = csv_input_file.columns.to_list()
    # if len(header) == 1:
    #     raise Exception("Invalid input file format: Input CSV file is not comma-separated.")
    def uniqueSNP(x):
        l = str(x['pileup'])
        # make upper
        l = l.upper()  # just to make simple
        l = utils.varOnly(l)

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

        df = pd.merge(so, df, on='snp_pos')
        return df


    # MAIN
    # mpileup input
    mp = pd.read_table(ARGS.in_mp, header=None)
    mp = mp.rename(columns={1: "snp_pos", 4: "pileup"})

    #snpeff-like table input
    so = pd.read_table(ARGS.in_so)

    sub_cols = ['chr_id', 'snp_pos', 'ref_allele',
                'alt_allele', 'codon1_genome_pos',
                'codon2_genome_pos', 'codon3_genome_pos',
                'ref_codon']

    so = so[sub_cols]

    mdf = utils.mergeTab(mp, so)

    df = mdf.apply(axis=1, func=uniqueSNP)
    df = pd.concat(list(df))
    reorg = ['chr_id', 'snp_pos', 'ref_allele', 'alt_allele',
             'varify_alt', 'codon1_genome_pos', 'codon2_genome_pos',
             'codon3_genome_pos', 'ref_codon', 'varify_codon1',
             'varify_codon2', 'varify_codon3', 'flag']

    df = df[reorg]

    # create report csv
    print(df)