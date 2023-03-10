# v.1.0.0
import argparse
import pandas as pd
import sys
import traceback
from pathlib import Path
import re

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

    # MAIN
    # mpileup input
    mp = pd.read_table(ARGS.in_mp, header=None)
    mp = mp.rename(columns={1: "snp_pos", 4: "pileup"})

    #snpeff-like table input
    so = pd.read_table(ARGS.in_so)

    utils.varify(mp, so)

    # create report csv