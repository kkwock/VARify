# NYU Bioinformatics Capstone Project: VARify
## Background and Purpose
With genetic testing becoming an integral component in clinical diagnosis, accuracy in its predictions are critical for proper assessment and treatment for patients. Current predictive tools use genomic data via whole exome sequencing (WES) and whole genome sequencing (WGS) to predict disease risk, however the diagnostic accuracy is ~40%, leaving more than half of patients without a proper diagnosis (Kim et al., 2020). While comprehensive, both methods lack functional evidence to properly classify disease-causing variants. 

For example, predicted effects of missense and nonsense mutations (i.e., premature termination codons) using tools like SnpEff (Cingolani et al., 2012) or VariantAnnotator (Obenchain et al., 2014) are based on gene model annotations from the reference genome. This can lead to inaccurate predictions in cases where a mutation is not expressed, or where the use of the reference annotation leads to an incorrect effect prediction (e.g., when a sample transcriptome utilizes alternate transcripts or splice-site junctions are shifted relative to the reference). A study has shown evidence of mutant transcripts escaping nonsense-mediated decay (NMD), with their translations not always leading to the expected loss of function (Anderson et al., 2017). It is therefore critical to be able to identify whether mutational effects predicted from WGS or WES are actually expressed in the transcriptome of interest.

To bridge this gap, researchers are incorporating RNA-seq as a way to detect transcriptional defects and increase diagnostic yield (Curry et al., 2021; Docking et al., 2021; Yépez et al., 2022). However, there is no one tool to assess this. One clinical study utilized multiple tools in order to compare the results of each method and improve patient stratification (Docking et al., 2021). This highlights the need for a method to validate predicted mutants in expressed transcripts. 

As an effort to fulfill this need, the purpose of our tool, VARify, is to compare and capture the predicted variants, based on WGS/WES data, in expressed transcripts. This would provide additional information that may increase mutant prediction accuracy and aid in classifying disease-causing variants. 

## Contents
* `report.qmd`
* `scripts/utility.R`
* `scripts/temp.sh`

# Tool Method Plan
![image](https://user-images.githubusercontent.com/92037458/210402025-93ef1bd0-7528-4f35-9754-23ad3988c401.png)

# References
* Anderson, J. L., Mulligan, T. S., Shen, M.-C., Wang, H., Scahill, C. M., Tan, F. J., Du, S. J., Busch-Nentwich, E. M., & Farber, S. A. (2017). mRNA processing in mutant zebrafish lines generated by chemical and CRISPR-mediated mutagenesis produces unexpected transcripts that escape nonsense-mediated decay. PLoS Genetics, 13(11), e1007105.
* Cingolani, P., Platts, A., Wang, L. L., Coon, M., Nguyen, T., Wang, L., Land, S. J., Lu, X., & Ruden, D. M. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly, 6(2), 80.
* Curry, P. D. K., Broda, K. L., & Carroll, C. J. (2021). The Role of RNA-Sequencing as a New Genetic Diagnosis Tool. Current Genetic Medicine Reports, 9(2), 13–21.
* Docking, T. R., Parker, J. D. K., Jädersten, M., Duns, G., Chang, L., Jiang, J., Pilsworth, J. A., Swanson, L. A., Chan, S. K., Chiu, R., Nip, K. M., Mar, S., Mo, A., Wang, X., Martinez-Høyer, S., Stubbins, R. J., Mungall, K. L., Mungall, A. J., Moore, R. A., … Karsan, A. (2021). A clinical transcriptome approach to patient stratification and therapy selection in acute myeloid leukemia. Nature Communications, 12(1), 2474.
* Kim, M.-J., Yum, M.-S., Seo, G. H., Lee, Y., Jang, H. N., Ko, T.-S., & Lee, B. H. (2020). Clinical Application of Whole Exome Sequencing to Identify Rare but Remediable Neurologic Disorders. Journal of Clinical Medicine Research, 9(11). https://doi.org/10.3390/jcm9113724
* Obenchain, V., Lawrence, M., Carey, V., Gogarten, S., Shannon, P., & Morgan, M. (2014). VariantAnnotation : a Bioconductor package for exploration and annotation of genetic variants. Bioinformatics , 30(14), 2076–2078.
* Yépez, V. A., Gusic, M., Kopajtich, R., Mertes, C., Smith, N. H., Alston, C. L., Ban, R., Beblo, S., Berutti, R., Blessing, H., Ciara, E., Distelmaier, F., Freisinger, P., Häberle, J., Hayflick, S. J., Hempel, M., Itkis, Y. S., Kishita, Y., Klopstock, T., … Prokisch, H. (2022). Clinical implementation of RNA sequencing for Mendelian disease diagnostics. Genome Medicine, 14(1), 38.
