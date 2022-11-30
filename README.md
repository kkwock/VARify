# NYU Bioinformatics Capstone Project
Repository for NYU Capstone Project with Dr. Flowers

# Background and Purpose
A nonsense mutation, or premature stop mutation, is a change in the DNA sequence that causes a protein to terminate earlier than expected. This would result in a truncated, and likely nonfunctional protein that can disrupt the normal regulation of an organism and cause disease. To mitigate this issue, eukaryotes have a nonsense-mediated decay (NMD) surveillance pathway that degrades mRNA transcripts with premature stop codons. 

However, a previous study has shown evidence of mutant transcripts escaping NMD and their translations did not always lead to loss of function (2017, Anderson et. al). This sheds light on the need for a predictive method for LOF with nonsense mutations. 

Current computational variant effect predictors, even those based solely on sequence conservation, underperform on non-LOF mutations (2022, Gerasimavicius et al.). The purpose of this tool is to identify early stop-codons and predict a loss-of-function for the protein. 

# Contents
* `report.qmd`
* `scripts/utility.R`
* `scripts/temp.sh`

# References
1. https://pubmed.ncbi.nlm.nih.gov/29161261/
2. https://www.nature.com/articles/s41598-020-71457-1
3. https://www.nature.com/articles/s41467-022-31686-6
