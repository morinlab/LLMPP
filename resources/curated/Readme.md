# Curated lymphoma gene lists

## Details

These lists are meant to be comprehensive and, as a result, contain many genes that are likely not relevant to these malignancies. Be sure to filter the lists based on how stringent you would like to be. Most people will be interested in the core lists. 

The master curated gene list for DLBCL including all genes nominated by any exome/genome-wide study can be found in this directory in `dlbcl_genes.tsv` or [here](dlbcl_genes.tsv). Most of the columns in this file are self-explanatory. The earliest_support column is meant to refer to the PubMed ID of the first study that nominated that gene as a significantly mutated gene in DLBCL. The columns Chapuy,	Reddy and	LymphGen indicate TRUE/FALSE for whether each gene was nominated/reported by that study. The next column (curated) is TRUE only for genes that have made it to the curated core list of DLBCL genes. The Lacy column indicates whether the gene was sequenced by Lacy et al. 

The master curated list for BL including DLBCL genes that have been scrutinized within BL is can be found in `bl_genes.tsv` (or [here](bl_genes.tsv)). The earliest_support_BL column indicates the PubMed ID of the first study to nominate this gene as mutated in BL (or NA when not applicable). Similar to dlbcl_genes.tsv, this file contains a Curated_BL_driver column, which will be TRUE for the subset of genes that have made it to our core list. The current core list at the time this document was created is below. This may not match the actual file, depending on whether this document is kept up to date. Please refer to the the file rather than this list. 

```
Gene	Curated_BL_driver	frequency_BL_Thomas
MYC	TRUE	60.2
DDX3X	TRUE	48.7
ID3	TRUE	47
TP53	TRUE	41.9
ARID1A	TRUE	36
CCND3	TRUE	28
FOXO1	TRUE	28
FBXO11	TRUE	21.6
GNA13	TRUE	21.6
SMARCA4	TRUE	17.8
KMT2D	TRUE	14
PCBP1	TRUE	12.3
P2RY8	TRUE	11
TCF3	TRUE	11
SIN3A	TRUE	10.6
TFAP4	TRUE	10.6
GNAI2	TRUE	9.7
RFX7	TRUE	9.3
BMP7	TRUE	8.9
CHD8	TRUE	8.5
RHOA	TRUE	8.1
EPPK1	TRUE	7.6
USP7	TRUE	6.8
WNK1	TRUE	6.8
HNRNPU	TRUE	6.4
PHF6	TRUE	5.5
EIF4A1	TRUE	5.1
PTEN	TRUE	4.7
```

## How to contribute

These lists may be incomplete and may contain errors. If you believe a gene is missing, is mis-attributed, etc please let us know. Feedback can be provided by submitting a GitHub issue. 

## How to cite

If you use the information in these lists please cite the following papers:

### DLBCL gene list

Minimum Information for Reporting a Genomics Experiment. Dreval K, Boutros PC, Morin RD.
Blood. 2022 Oct 11:blood.2022016095. doi: 10.1182/blood.2022016095. PMID: 36219881

### BL gene list

GENETIC SUBGROUPS INFORM ON PATHOBIOLOGY IN ADULT AND PEDIATRIC BURKITT LYMPHOMA. Thomas N, Dreval K et al. Blood. 2022 Oct 6:blood.2022016534. doi: 10.1182/blood.2022016534. PMID: 36201743
