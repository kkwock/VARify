import pysam

# Replace "input.bam" with the actual path to your BAM file
bamfile = "/Users/kmae/Documents/GitHub/VARify/data/test/ENCFF283TLK_chr12_pos40509062.txt"

# Open the BAM file
with pysam.AlignmentFile(bamfile, "rb") as samfile:
    # Iterate through alignments
    for alignment in samfile:
        # Get desired fields from the alignment
        chromosome = alignment.reference_name
        position = alignment.reference_start + 1  # Adding 1 to convert to 1-based position
        md_tag = alignment.get_tag("MD", with_value_type=True)
        cigar_string = alignment.cigarstring

        # Print or process the extracted information
        print(f"Chromosome: {chromosome}, Position: {position}, MD Tag: {md_tag}, CIGAR: {cigar_string}")