# Curated lymphoma gene lists

# What these files are and are not

These lists are meant to be comprehensive lists of genes that are (reportedly) enriched for non-silent or other driver mutations in one or more B-cell lymphomas. As a result, contain many genes that are likely not relevant to these malignancies and some genes may be missing if they are only affected by aberrant somatic hypermutation (aSHM) but not demonstrated to be _bona fide_ drivers. Comprehensive lists of regions affected by aSHM are provided separately but are not considered gene lists per-se because they are not all within genes. The completeness of these lists is a blessing and a curse. Be sure to filter the lists based on how stringent you would like to be. Most people will be interested in just the core lists.

## The lists

The master curated gene list for DLBCL including all genes nominated by any exome/genome-wide study can be found in this directory in `dlbcl_genes.tsv` or [here](dlbcl_genes.tsv). Most of the columns in this file are self-explanatory. The earliest_support column is meant to refer to the PubMed ID of the first study that nominated that gene as a significantly mutated gene in DLBCL. The columns Chapuy,	Reddy and	LymphGen indicate TRUE/FALSE for whether each gene was nominated/reported by that study. The next column (curated) is TRUE only for genes that have made it to the curated core list of DLBCL genes. The Lacy column indicates whether the gene was sequenced by Lacy et al. The core list at the time this document was prepared is shown below along with a few key columns from the file. _This may not match the actual file, depending on whether this document is kept up to date. Please refer to the the file rather than this list._

| Gene      | curated | aSHM  | earliest_support |
|-----------|---------|-------|------------------|
| ACTB      | TRUE    | TRUE  | 22343534         |
| ARID1A    | TRUE    | FALSE | 23292937         |
| ATM       | TRUE    | FALSE | 28985567         |
| B2M       | TRUE    | FALSE | 21796119         |
| BCL10     | TRUE    | FALSE |                  |
| BCL11A    | TRUE    | TRUE  |                  |
| BCL2      | TRUE    | TRUE  |                  |
| BCL6      | TRUE    | TRUE  | 21796119         |
| BCL7A     | TRUE    | TRUE  |                  |
| BCOR      | TRUE    | FALSE |                  |
| BCR       | TRUE    | FALSE |                  |
| BIRC6     | TRUE    | FALSE | 28985567         |
| BRAF      | TRUE    | FALSE | 22343534         |
| BTG1      | TRUE    | TRUE  | 21796119         |
| BTG2      | TRUE    | TRUE  | 21796119         |
| BTK       | TRUE    | FALSE |                  |
| CARD11    | TRUE    | FALSE | 21796119         |
| CASP8     | TRUE    | FALSE | 28985567         |
| CCND3     | TRUE    | FALSE | 21796119         |
| CD274     | TRUE    | FALSE | 21796119         |
| CD36      | TRUE    | FALSE |                  |
| CD58      | TRUE    | FALSE | 21796119         |
| CD70      | TRUE    | FALSE | 21796119         |
| CD79B     | TRUE    | FALSE | 21796119         |
| CD83      | TRUE    | TRUE  | 29641966         |
| CDKN2A    | TRUE    | FALSE |                  |
| CIITA     | TRUE    | TRUE  | 21796119         |
| CREBBP    | TRUE    | FALSE |                  |
| CXCR4     | TRUE    | TRUE  | 23131835         |
| DAZAP1    | TRUE    | FALSE |                  |
| DDX3X     | TRUE    | FALSE | 29641966         |
| DTX1      | TRUE    | TRUE  | 29641966         |
| DUSP2     | TRUE    | TRUE  | 28985567         |
| EBF1      | TRUE    | TRUE  | 23174882         |
| EP300     | TRUE    | FALSE |                  |
| ETS1      | TRUE    | TRUE  | 21796119         |
| ETV6      | TRUE    | TRUE  |                  |
| EZH2      | TRUE    | FALSE | 20081860         |
| FAS       | TRUE    | FALSE |                  |
| FBXO11    | TRUE    | FALSE |                  |
| FBXW7     | TRUE    | FALSE |                  |
| FOXC1     | TRUE    | FALSE |                  |
| FOXO1     | TRUE    | FALSE | 21796119         |
| GNA13     | TRUE    | FALSE | 21796119         |
| GNAI2     | TRUE    | FALSE |                  |
| GRHPR     | TRUE    | TRUE  |                  |
| HIST1H1B  | TRUE    | TRUE  |                  |
| HIST1H1C  | TRUE    | TRUE  | 21796119         |
| HIST1H1D  | TRUE    | TRUE  |                  |
| HIST1H1E  | TRUE    | TRUE  |                  |
| HIST1H2AC | TRUE    | TRUE  |                  |
| HIST1H2AM | TRUE    | TRUE  |                  |
| HIST1H2BC | TRUE    | TRUE  | 28985567         |
| HIST1H2BK | TRUE    | TRUE  |                  |
| HIST1H3B  | TRUE    | TRUE  |                  |
| HIST2H2BE | TRUE    | TRUE  |                  |
| HLA-A     | TRUE    | FALSE |                  |
| HLA-B     | TRUE    | FALSE |                  |
| HLA-C     | TRUE    | FALSE |                  |
| HLA-DMA   | TRUE    | FALSE |                  |
| HLA-DMB   | TRUE    | FALSE |                  |
| HNRNPD    | TRUE    | FALSE |                  |
| HNRNPH1   | TRUE    | FALSE |                  |
| HNRNPU    | TRUE    | FALSE |                  |
| HVCN1     | TRUE    | FALSE |                  |
| ID3       | TRUE    | FALSE | 22885699         |
| IL16      | TRUE    | FALSE |                  |
| IL4R      | TRUE    | TRUE  | 33684939         |
| IRF4      | TRUE    | TRUE  | 21796119         |
| IRF8      | TRUE    | TRUE  | 21796119         |
| ITPKB     | TRUE    | TRUE  | 29641966         |
| KLF2      | TRUE    | TRUE  | 29641966         |
| KLHL14    | TRUE    | FALSE | 23292937         |
| KLHL6     | TRUE    | TRUE  | 21796119         |
| KMT2C     | TRUE    | FALSE | 23292937         |
| KMT2D     | TRUE    | FALSE | 21796119         |
| KRAS      | TRUE    | FALSE | 22343534         |
| LCOR      | TRUE    | FALSE |                  |
| LTB       | TRUE    | TRUE  |                  |
| LYN       | TRUE    | FALSE |                  |
| MCL1      | TRUE    | FALSE | 28985567         |
| MEF2B     | TRUE    | TRUE  | 21796119         |
| MEF2C     | TRUE    | TRUE  |                  |
| MGA       | TRUE    | FALSE | 23292937         |
| MPEG1     | TRUE    | FALSE |                  |
| MS4A1     | TRUE    | TRUE  |                  |
| MTOR      | TRUE    | FALSE | 23292937         |
| MYC       | TRUE    | TRUE  |                  |
| MYD88     | TRUE    | FALSE | 21179087         |
| MYOM2     | TRUE    | FALSE |                  |
| NCOR2     | TRUE    | FALSE |                  |
| NFKB1     | TRUE    | FALSE |                  |
| NFKBIA    | TRUE    | FALSE |                  |
| NFKBIE    | TRUE    | FALSE |                  |
| NFKBIZ    | TRUE    | TRUE  |                  |
| NLRC5     | TRUE    | FALSE |                  |
| NOL9      | TRUE    | FALSE |                  |
| NOTCH1    | TRUE    | FALSE | 22343534         |
| NOTCH2    | TRUE    | FALSE |                  |
| OSBPL10   | TRUE    | TRUE  |                  |
| P2RY8     | TRUE    | FALSE | 22343534         |
| PCBP1     | TRUE    | FALSE |                  |
| PIM1      | TRUE    | TRUE  | 11460166         |
| PIM2      | TRUE    | TRUE  |                  |
| POU2AF1   | TRUE    | TRUE  |                  |
| POU2F2    | TRUE    | FALSE |                  |
| PPP1R9B   | TRUE    | FALSE |                  |
| PRDM1     | TRUE    | FALSE |                  |
| PRKDC     | TRUE    | FALSE |                  |
| PTEN      | TRUE    | FALSE |                  |
| PTPRD     | TRUE    | FALSE |                  |
| RB1       | TRUE    | FALSE |                  |
| RFX7      | TRUE    | FALSE |                  |
| RFXAP     | TRUE    | FALSE |                  |
| RHOA      | TRUE    | FALSE | 11460166         |
| RRAGC     | TRUE    | FALSE | 26691987         |
| S1PR2     | TRUE    | TRUE  |                  |
| SETD1B    | TRUE    | FALSE |                  |
| SETD2     | TRUE    | FALSE |                  |
| SF3B1     | TRUE    | FALSE |                  |
| SGK1      | TRUE    | TRUE  | 21796119         |
| SIN3A     | TRUE    | FALSE |                  |
| SMARCA4   | TRUE    | FALSE | 23292937         |
| SOCS1     | TRUE    | TRUE  |                  |
| SPEN      | TRUE    | FALSE |                  |
| STAT3     | TRUE    | FALSE |                  |
| STAT6     | TRUE    | FALSE |                  |
| TBL1XR1   | TRUE    | FALSE |                  |
| TCL1A     | TRUE    | TRUE  |                  |
| TET2      | TRUE    | FALSE |                  |
| TMEM30A   | TRUE    | FALSE | 21796119         |
| TMSB4X    | TRUE    | TRUE  |                  |
| TNFAIP3   | TRUE    | FALSE |                  |
| TNFRSF14  | TRUE    | FALSE | 20884631         |
| TOX       | TRUE    | FALSE |                  |
| TP53      | TRUE    | FALSE |                  |
| TRAF3     | TRUE    | FALSE |                  |
| TRIP12    | TRUE    | FALSE |                  |
| TRRAP     | TRUE    | FALSE |                  |
| UBE2A     | TRUE    | FALSE |                  |
| UNC5C     | TRUE    | FALSE |                  |
| UNC5D     | TRUE    | FALSE |                  |
| USP7      | TRUE    | FALSE |                  |
| VPS13B    | TRUE    | FALSE |                  |
| WEE1      | TRUE    | FALSE |                  |
| XBP1      | TRUE    | FALSE |                  |
| XPO1      | TRUE    | FALSE | 26608593         |
| ZC3H12A   | TRUE    | FALSE |                  |
| ZFP36L1   | TRUE    | TRUE  |                  |
| ZNF292    | TRUE    | FALSE | 23292937         |
| CXCR5     | TRUE    | FALSE | 29641966         |
| TAP1      | TRUE    | FALSE | 29641966         |

The master curated list for BL including DLBCL genes that have been scrutinized within BL is can be found in `bl_genes.tsv` (or [here](bl_genes.tsv)). The earliest_support_BL column indicates the PubMed ID of the first study to nominate this gene as mutated in BL (or NA when not applicable). The frequency_BL_Thomas and frequency_BL_Panea columns report the percentage of patient samples with at least one non-silent mutation in this gene in the two studies. In the case of Panea et al, these numbers are based on the reanalysis of the exome data from this study by the Morin lab (as detailed in Dreval et al). The original frequencies based on the mutation calls from Panea et al are in the frequency_BL_Panea_original column.

Similar to dlbcl_genes.tsv, this file contains a Curated_BL_driver column, which will be TRUE for the subset of genes that have made it to our core list. The current core list at the time this document was created is below. _This may not match the actual file, depending on whether this document is kept up to date. Please refer to the the file rather than this list._ 

| Gene    | Curated_BL_driver | frequency_BL_Thomas | frequency_BL_Panea |
|---------|-------------------|---------------------|--------------------|
| MYC     | TRUE              | 60.2                | 49.5               |
| DDX3X   | TRUE              | 48.7                | 39.6               |
| ID3     | TRUE              | 47                  | 31.7               |
| TP53    | TRUE              | 41.9                | 43.6               |
| ARID1A  | TRUE              | 36                  | 19.8               |
| CCND3   | TRUE              | 28                  | 17.8               |
| FOXO1   | TRUE              | 28                  | 21.8               |
| FBXO11  | TRUE              | 21.6                | 15.8               |
| GNA13   | TRUE              | 21.6                | 20.8               |
| SMARCA4 | TRUE              | 17.8                | 18.8               |
| KMT2D   | TRUE              | 14                  | 15.8               |
| PCBP1   | TRUE              | 12.3                | 11.9               |
| P2RY8   | TRUE              | 11                  | NA                 |
| TCF3    | TRUE              | 11                  | 9.9                |
| SIN3A   | TRUE              | 10.6                | 14.9               |
| TFAP4   | TRUE              | 10.6                | 9.9                |
| GNAI2   | TRUE              | 9.7                 | 8.9                |
| RFX7    | TRUE              | 9.3                 | 4                  |
| BMP7    | TRUE              | 8.9                 | 5                  |
| CHD8    | TRUE              | 8.5                 | 16.8               |
| RHOA    | TRUE              | 8.1                 | 12.9               |
| EPPK1   | TRUE              | 7.6                 | 15.8               |
| USP7    | TRUE              | 6.8                 | 5.9                |
| WNK1    | TRUE              | 6.8                 | 11.9               |
| HNRNPU  | TRUE              | 6.4                 | 8.9                |
| PHF6    | TRUE              | 5.5                 | 5.9                |
| EIF4A1  | TRUE              | 5.1                 | 6.9                |
| PTEN    | TRUE              | 4.7                 | 4                  |

## How to contribute

These lists may be incomplete and may contain errors. If you believe a gene is missing, is mis-attributed, etc please let us know. Feedback can be provided by submitting a GitHub issue. 

## How to cite

If you use the information in these lists please cite the following papers:

### DLBCL gene list

Minimum Information for Reporting a Genomics Experiment. Dreval K, Boutros PC, Morin RD.
Blood. 2022 Oct 11:blood.2022016095. doi: 10.1182/blood.2022016095. PMID: 36219881

### BL gene list

GENETIC SUBGROUPS INFORM ON PATHOBIOLOGY IN ADULT AND PEDIATRIC BURKITT LYMPHOMA. Thomas N, Dreval K et al. Blood. 2022 Oct 6:blood.2022016534. doi: 10.1182/blood.2022016534. PMID: 36201743
