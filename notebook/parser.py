import argparse
import pandas as pd
from multiprocessing import Pool, cpu_count

# Define the remove_softclip function
def remove_softclip(s):
    while s.startswith('N'):
        s = s[1:]
    return s

def process_file(input_file):
    # Read the CSV file into a DataFrame
    s2p_out = pd.read_csv(input_file, header=None)

    # Extract every 1st and 4th entry into separate DataFrames
    col1 = s2p_out.iloc[::4].reset_index(drop=True)
    col2 = s2p_out.iloc[3::4].reset_index(drop=True)

    # Split the data in column 1 by tab delimiter and keep elements at indices 1, 3, 4, and 9
    split_col1 = col1[0].str.split('\t', expand=True).iloc[:, [0, 2, 3, 8]]

    # Apply remove_softclip function to each record in col2
    col2_cleaned = col2.apply(lambda x: remove_softclip(x[0]), axis=1)

    # Create a new DataFrame with the split data and col2_cleaned as columns
    result_df = pd.concat([split_col1, col2_cleaned], axis=1)

    # Rename the columns
    result_df.columns = ['QNAME', 'RNAME', 'POS', 'TLEN', 'SEQ']

    # Return the result DataFrame
    return result_df

def main(input_file):
    # Determine the number of CPUs to use for multiprocessing
    num_cpus = cpu_count()

    # Create a Pool of workers
    pool = Pool(processes=num_cpus)

    # Process the input file using multiple processes
    result_dfs = pool.map(process_file, [input_file])

    # Close the pool of workers
    pool.close()

    # Concatenate the result DataFrames (there's only one in this case)
    final_result = pd.concat(result_dfs)

    # Print the final result DataFrame
    return final_result

if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description="Process sam2pair_5Entry.out file")

    # Add arguments
    parser.add_argument("input_file", type=str, help="Path to the input *.out file")

    # Parse arguments
    args = parser.parse_args()

    # Call the main function with input_file as argument
    main(args.input_file)
