# AtDermQTLGWAS_manuscript

Repository to centralize code, figures and material related to the "Integration of GWAS, QTLs and keratinocyte functional assays reveals molecular mechanisms of Atopic Dermatitis" Oliva et al 2025 . It is structured as follows:

# 1. Introduction
Atopic dermatitis is a highly heritable and common inflammatory skin condition affecting children and adults worldwide. Multi-ancestry approaches to atopic dermatitis genetic association studies are poised to boost power to detect genetic signal and identify loci contributing to atopic dermatitis risk. Here, we present a multi-ancestry GWAS meta-analysis of twelve atopic dermatitis cohorts from five ancestral populations totaling 56,146 cases and 602,280 controls. We report 101 genomic loci associated with atopic dermatitis, including 16 loci that have not been previously associated with atopic dermatitis or eczema. Fine-mapping, QTL colocalization, and cell-type enrichment analyses identified genes and cell types implicated in atopic dermatitis pathophysiology. Functional analyses in keratinocytes provide evidence for genes that could play a role in atopic dermatitis through epidermal barrier function. Our study provides insights into the etiology of atopic dermatitis by harnessing multiple genetic and functional approaches to unveil the mechanisms by which atopic dermatitis-associated variants impact genes and cell types.

# 2. Contents

## Analyses
. One folder per main analysis.

. Each folder needs to include a) ioptionally, a file named 'README', with instructions of how to execute the code, b) a folder named 'code', containing i) one (or more) master script files or juniper notebook files with code to run analyses sequentially and ii) scripts containing the code, called in (i), and c) a folder named 'data'. This folder needs to contain references to the data that serves as input to (b) scripts.

e.g.
### cell_type_enrichment_MAGMA
Includes code to perform GWAS SNP enrichment analysis to connect GWAS signals with cell-type specific functional annotations and hence, priritize the cell types of origin of AD.

### GWAS_QTL_colocalization
Includes code for QTL-GWAS colocalization analysis using coloc. coloc is available at https://github.com/chr1swallace/coloc/.

## Figures
Includes code and data to generate manuscript figures. Additional data can be found at https://doi.org/10.25452/figshare.plus.c.7282969

. One folder per figure.

. Each folder needs to include a) a folder named 'code', containing scripts with the code to generate the figure, and b) a folder named 'data',  containing processed data to generate the figure. If data is too big to be included in the github, the 'data' folder needs to contain symbolic links / paths / sources to the data that serves as input to the figure.
