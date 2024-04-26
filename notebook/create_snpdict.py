'''
create_snpdict.py
description: creates dictionary of codons
input: sam2pairwise output table, snpeff-like table
output: codons.json
'''

# read in data
s2p = pd.read_csv(s2pout, header=None) # need to ensure None header to work
snpeff_table = pd.read_csv(snpeff, sep='\t')

# extract sam2pairwise data
result_df = sam2pair_extract(s2p)

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

# Writing to codons.json - provides raw codons
with open("codons.json", "w") as outfile:
    outfile.write(codons)
