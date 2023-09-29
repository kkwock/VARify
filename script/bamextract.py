import pysam
import re

# Replace "input.bam" with the actual path to your BAM file
bamfile = "/Users/kmae/Documents/GitHub/VARify/data/test/ENCFF283TLK_chr12_pos40509062.txt"

# # Open the BAM file
# with pysam.AlignmentFile(bamfile, "rb") as samfile:
#     # Iterate through alignments
#     for alignment in samfile:
#         # Get desired fields from the alignment
#         chromosome = alignment.reference_name
#         position = alignment.reference_start + 1  # Adding 1 to convert to 1-based position
#         md_tag = alignment.get_tag("MD", with_value_type=True)
#         cigar_string = alignment.cigarstring
#
#         # Print or process the extracted information
#         print(f"Chromosome: {chromosome}, Position: {position}, MD Tag: {md_tag}, CIGAR: {cigar_string}")


# Open the SAM-formatted text file
with open(bamfile, "r") as file:
    content = file.read()

    # Split the content into a list of strings
    fields = re.split(r'\s+', content)

    # Now you can access individual fields using indexing
    print(fields[0])  # Print the first field
    print(fields[1])  # Print the second field

with open(bamfile, "r") as file:
    # Split the content into a list of strings
    content = file.read()
    fields = re.split(r'\s+', content)

    for line in file:
        if line.startswith("@"):
            # Skip header lines
            continue

        fields = line.strip().split("\t")
        print(fields)  # Print the split fields to diagnose the issue
        chromosome = fields[2]
        position = int(fields[3])

        tags = fields[11:]
        md_tag = None
        for tag in tags:
            if tag.startswith("MD:Z:"):
                md_tag = tag.split(":")[2]
                break

        cigar_string = fields[5]

        print(f"Chromosome: {chromosome}, Position: {position}, MD Tag: {md_tag}, CIGAR: {cigar_string}")