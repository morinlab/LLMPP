# Curated lymphoma gene lists

# What these files are and what they are not

These lists are meant to be comprehensive lists of genes that are (reportedly) enriched for non-silent or other driver mutations in one or more B-cell lymphomas. As a result, contain many genes that are likely not relevant to these malignancies and some genes may be missing if they are only affected by aberrant somatic hypermutation (aSHM) but not demonstrated to be _bona fide_ drivers. Comprehensive lists of regions affected by aSHM are provided separately but are not considered gene lists per-se because they are not all within genes. The completeness of these lists is a blessing and a curse. Be sure to filter the lists based on how stringent you would like to be. Most people will be interested in just the core lists.

## The lists

### DLBCL

The master curated gene list for DLBCL including all genes nominated by any exome/genome-wide study can be found in this directory in `dlbcl_genes.tsv` or [here](dlbcl_genes.tsv). Most of the columns in this file are self-explanatory. The earliest_support column is meant to refer to the PubMed ID of the first study that nominated that gene as a significantly mutated gene in DLBCL. The columns Chapuy,	Reddy and	LymphGen indicate TRUE/FALSE for whether each gene was nominated/reported by that study. The next column (curated) is TRUE only for genes that have made it to the curated core list of DLBCL genes. 
This file contains four additional columns, which respectively report the mutation frequency in the cohorts from Chapuy et al, Schmitz et al, Reddy et al If a gene has been identified from additional information, this is from our compendium of genomes in the GAMBL project. The Lacy column indicates whether the gene was sequenced by Lacy et al. The core gene list at the time this document was prepared is shown below along with a few key columns from the file. _This may not match the actual file, depending on whether this document is kept up to date. Please refer to the the file rather than this list._ The separate file `dlbcl_genes_with_mutation_frequencies.tsv` ([here](dlbcl_genes_with_mutation_frequencies.tsv)) is intended to include all genes nominated to be recurrently/significantly mutated in DLBCL. It contains columns indicating, where possible, which study nominated the gene. 



|Gene|Tier|aSHM|QC|Mean variant quality|citekey|PMID|MutationEffect|Mutation-PMID|MutationEffect-citekey|LymphGen|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|ACTB|1|TRUE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534||||TRUE
|ACTG1|1|TRUE|NA||fanComprehensiveCharacterizationDriver2020|32565964||||FALSE
|ARID1A|1|FALSE|NA||zhangGeneticHeterogeneityDiffuse2013|23292937|LOF|38458187|barisicARID1AOrchestratesSWI2024|TRUE
|ATM|1|FALSE|NA||reddyGeneticFunctionalDrivers2017|28985567|LOF|11756177|camachoATMGeneInactivation2002|FALSE
|B2M|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|22137796|challa-malladiCombinedGeneticInactivation2011|FALSE
|BCL10|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|GOF|35658124|xiaBCL10MutationsDefine2022|TRUE
|BCL2|1|TRUE|NA||tanakaFrequentIncidenceSomatic1992|1339299||||TRUE
|BCL6|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|12504096|masclePointMutationsBCL62003|TRUE
|BCL7A|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|32576963|balinas-gaviraFrequentMutationsAminoterminal2020|FALSE
|BIRC6|1|FALSE|NA||reddyGeneticFunctionalDrivers2017|28985567||||FALSE
|BRAF|1|FALSE|NA||tiacciBRAFMutationsHairycell2011|22343534|LOF|15035987|wanMechanismActivationRAFERK2004|FALSE
|BTG1|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|33021411|almasmoumFrequentLossBTG12021|TRUE
|BTG2|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119||||TRUE
|BTK|1|FALSE|NA||albuquerqueEnhancingKnowledgeDiscovery2017|28327945|LOF|33419778|huFollicularLymphomaassociatedBTK2021|FALSE
|CARD11|1|FALSE|NA||lenzOncogenicCARD11Mutations2008|18323416|GOF|18323416|lenzOncogenicCARD11Mutations2008|FALSE
|CCND3|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|GOF|22885699|schmitzBurkittLymphomaPathogenesis2012|FALSE
|CD58|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|22137796|challa-malladiCombinedGeneticInactivation2011|TRUE
|CD70|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|36471481|nieDualRoleCD702022|TRUE
|CD79B|1|FALSE|NA||davisChronicActiveBcellreceptor2010|20054396|GOF|20054396|davisChronicActiveBcellreceptor2010|TRUE
|CD83|1|TRUE|NA||morinMutationalStructuralAnalysis2013|23699601||||TRUE
|CDKN2A|1|FALSE|NA||morinMutationalStructuralAnalysis2013|23699601|LOF|19260062|kannengiesserFunctionalStructuralGenetic2009|TRUE
|CIITA|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|26549456|mottokGenomicAlterationsCIITA2015|TRUE
|CREBBP|1|FALSE|NA||pasqualucciInactivatingMutationsAcetyltransferase2011|21390126|LOF|21390126|pasqualucciInactivatingMutationsAcetyltransferase2011|TRUE
|CXCR4|1|TRUE|NA||khodabakhshiRecurrentTargetsAberrant2012|23131835|LOF|36089616|zmajkovicovaGenotypephenotypeCorrelationsWHIM2022|FALSE
|DDX3X|1|FALSE|NA||reddyGeneticFunctionalDrivers2017|28985567|LOF|34437837|gongSequentialInverseDysregulation2021|TRUE
|DTX1|1|TRUE|NA||zhangGeneticHeterogeneityDiffuse2013|23292937||||TRUE
|DUSP2|1|TRUE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534||||TRUE
|EBF1|1|TRUE|NA||bohleRoleEarlyBcell2013|23174882|LOF|28692033|ramirez-komoSpontaneousLossLineage2017|FALSE
|EEF1A1|1|FALSE|NA||chapuyMolecularSubtypesDiffuse2018|29713087||||FALSE
|EP300|1|FALSE|NA||pasqualucciInactivatingMutationsAcetyltransferase2011|21390126|LOF|21390126|pasqualucciInactivatingMutationsAcetyltransferase2011|TRUE
|ETS1|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119||||TRUE
|ETV6|1|TRUE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534|LOF|24997145|wangETV6MutationCohort2014|TRUE
|EZH2|1|FALSE|NA||morinSomaticMutationsAltering2010|20081860|GOF|21078963|sneeringerCoordinatedActivitiesWildtype2010|TRUE
|FAS|1|FALSE|NA||schollMutationsRegionFAS2007|17487740|LOF|20935634|wangFasFADDDeathDomain2010|FALSE
|FBXO11|1|FALSE|NA||arthurGenomewideDiscoverySomatic2018|30275490|LOF|22113614|duanFBXO11TargetsBCL62011|FALSE
|FBXW7|1|FALSE|NA||zhangGeneticHeterogeneityDiffuse2013|23292937|LOF|32350066|saffieFBXW7TriggersDegradation2020|FALSE
|FOXO1|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|GOF|23460611|trinhAnalysisFOXO1Mutations|FALSE
|GNA13|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|25274307|muppidiLossSignalingGa132014|FALSE
|GNAI2|1|FALSE|NA||morinMutationalStructuralAnalysis2013|23699601|GOF|25274307|muppidiLossSignalingGa132014|FALSE
|GRB2|1|FALSE|NA||pasqualucciAnalysisCodingGenome2011|21804550||||FALSE
|GRHPR|1|TRUE|NA||schmitzGeneticsPathogenesisDiffuse2018|29641966||||TRUE
|HIST1H1B|1|TRUE|NA||chapuyMolecularSubtypesDiffuse2018|29713087||||FALSE
|HIST1H1C|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119||||FALSE
|HIST1H1D|1|TRUE|NA||morinMutationalStructuralAnalysis2013|23699601||||FALSE
|HIST1H1E|1|TRUE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534||||FALSE
|HIST1H2AC|1|TRUE|NA||morinMutationalStructuralAnalysis2013|23699601||||FALSE
|HIST1H2AM|1|TRUE|NA||chapuyMolecularSubtypesDiffuse2018|29713087||||FALSE
|HIST1H2BC|1|TRUE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534||||FALSE
|HIST1H2BK|1|TRUE|NA||zhangGeneticHeterogeneityDiffuse2013|23292937||||FALSE
|HIST1H3B|1|TRUE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534||||FALSE
|HIST2H2BE|1|TRUE|NA||chapuyMolecularSubtypesDiffuse2018|29713087||||FALSE
|HLA-A|1|FALSE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534|LOF|34050029|fangazioGeneticMechanismsHLAI2021|TRUE
|HLA-B|1|FALSE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534|LOF|34050029|fangazioGeneticMechanismsHLAI2021|TRUE
|HLA-C|1|FALSE|NA||chapuyMolecularSubtypesDiffuse2018|29713087|LOF|34050029|fangazioGeneticMechanismsHLAI2021|FALSE
|HLA-DMB|1|FALSE|NA||hubschmannMutationalMechanismsShaping2021|33953289||||FALSE
|HNRNPU|1|FALSE|NA||reddyGeneticFunctionalDrivers2017|28985567||||FALSE
|HVCN1|1|FALSE|NA||chapuyMolecularSubtypesDiffuse2018|29713087||||FALSE
|IKZF3|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|GOF|33689703|lazarianHotspotMutationTranscription2021|FALSE
|IL4R|1|TRUE|NA||dunsCharacterizationDLBCLPMBL2021|33684939|GOF|29467182|viganoSomaticIL4RMutations2018|FALSE
|IRF4|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119||||TRUE
|IRF8|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|38996030|qiuIRF8mutantCellLymphoma2024|TRUE
|ITPKB|1|TRUE|NA||schmitzGeneticsPathogenesisDiffuse2018|29641966|LOF|29650799|tiacciPervasiveMutationsJAKSTAT2018|TRUE
|JUNB|1|FALSE|PASS|3|lohrDiscoveryPrioritizationSomatic2012|22343534||||TRUE
|KLF2|1|TRUE|NA||pasqualucciAnalysisCodingGenome2011|21804550||||TRUE
|KLHL14|1|FALSE|NA||zhangGeneticHeterogeneityDiffuse2013|23292937|LOF|32127472|choiRegulationCellReceptordependent2020|TRUE
|KLHL6|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|29695787|choiLossKLHL6Promotes2018|FALSE
|KMT2C|1|FALSE|NA||zhangGeneticHeterogeneityDiffuse2013|23292937||||FALSE
|KMT2D|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|26366712|zhangDisruptionKMT2DPerturbs2015|TRUE
|KRAS|1|FALSE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534|GOF|9219684|scheffzekRasRasGAPComplexStructural1997|FALSE
|LRRN3|1|FALSE|NA||zhangGeneticHeterogeneityDiffuse2013|23292937||||FALSE
|LTB|1|TRUE|NA||chapuyMolecularSubtypesDiffuse2018|29713087||||FALSE
|MEF2B|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119|NEO|23974956; 26245647|ponMEF2BMutationsNonHodgkin2015|TRUE
|MEF2C|1|TRUE|NA||hubschmannMutationalMechanismsShaping2021|33953289||||FALSE
|MGA|1|FALSE|NA||reddyGeneticFunctionalDrivers2017|28985567|LOF|23039309|depaoliMGASuppressorMYC2013|FALSE
|MIR142|1|TRUE|NA||kwanhianMicroRNA142Mutated202012|23342264|LOF|29724719|trissalMIR142LossofFunctionMutations2018|FALSE
|MPEG1|1|FALSE|NA||morinMutationalStructuralAnalysis2013|23699601||||TRUE
|MS4A1|1|TRUE|NA||rushtonGeneticEvolutionaryPatterns2020|32589730|LOF|32589730|rushtonGeneticEvolutionaryPatterns2020|FALSE
|MTOR|1|FALSE|NA||zhangGeneticHeterogeneityDiffuse2013|23292937|GOF|24631838|grabinerDiverseArrayCancerassociated2014|FALSE
|MYC|1|TRUE|NA||pasqualucciHypermutationMultipleProtooncogenes2001|11460166|GOF|38565249|freieGermlinePointMutation2024|FALSE
|MYD88|1|FALSE|NA||ngoOncogenicallyActiveMYD882011|21179087|GOF|21179087|ngoOncogenicallyActiveMYD882011|TRUE
|NFKBIA|1|FALSE|NA||thomasMutationalAnalysisIkappaBalpha2004|15198731|LOF|10637284|jungnickelClonalDeleteriousMutations2000|TRUE
|NFKBIE|1|FALSE|NA||morinGeneticLandscapesRelapsed2016|26647218|LOF|25987724|mansouriFunctionalLossIkBe2015|FALSE
|NFKBIZ|1|TRUE|NA||morinGeneticLandscapesRelapsed2016|26647218|GOF|302754900|arthurGenomewideDiscoverySomatic2018|FALSE
|NOL9|1|FALSE|NA||schmitzGeneticsPathogenesisDiffuse2018|29641966||||TRUE
|NOTCH1|1|FALSE|NA||pasqualucciAnalysisCodingGenome2011|21804550|GOF|29045844|ryanCellRegulomeLinks2017|TRUE
|NOTCH2|1|FALSE|NA||zhangGeneticHeterogeneityDiffuse2013|23292937|GOF|19445024|leeGainoffunctionMutationsCopy2009|TRUE
|OSBPL10|1|TRUE|NA||zhangGeneticHeterogeneityDiffuse2013|23292937||||TRUE
|P2RY8|1|FALSE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534|LOF|25274307|muppidiLossSignalingGa132014|FALSE
|PIM1|1|TRUE|NA||pasqualucciHypermutationMultipleProtooncogenes2001|11460166|GOF|27904766|kuoRolePIM1Ibrutinibresistant2016|TRUE
|PIM2|1|TRUE|NA||reddyGeneticFunctionalDrivers2017|28985567||||TRUE
|POU2AF1|1|TRUE|NA||chapuyMolecularSubtypesDiffuse2018|29713087|LOF|30802265|gonzalez-rinconUnravelingTransformationFollicular2019|FALSE
|POU2F2|1|FALSE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534|LOF|26993806|hodsonRegulationNormalBcell2016|FALSE
|PRDM1|1|FALSE|NA||pasqualucciInactivationPRDM1BLIMP12006|NA|LOF|16492805|pasqualucciInactivationPRDM1BLIMP12006|TRUE
|PTEN|1|FALSE|NA||reddyGeneticFunctionalDrivers2017|28985567|LOF|23840064|pfeiferPTENLossDefines2013|FALSE
|PTPN6|1|FALSE|PASS|3|reddyGeneticFunctionalDrivers2017|28985567|LOF|26565811|demosthenousLossFunctionMutations2015|FALSE
|RB1|1|FALSE|NA||morinMutationalStructuralAnalysis2013|23699601|LOF|17332242|pinyolInactivationRB1Mantlecell2007|FALSE
|RFX7|1|FALSE|NA||arthurGenomewideDiscoverySomatic2018|30275490|LOF|30926791|weberPiggyBacTransposonTools2019|FALSE
|RHOA|1|FALSE|NA||zhangGeneticHeterogeneityDiffuse2013|23292937|LOF|26616858|ohayreInactivatingMutationsGNA132016|FALSE
|RRAGC|1|FALSE|NA||okosunRecurrentMTORC1activatingRRAGC2016|26691987|GOF|26691987|okosunRecurrentMTORC1activatingRRAGC2016|FALSE
|S1PR2|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|25274307|muppidiLossSignalingGa132014|TRUE
|SETD1B|1|FALSE|NA||albuquerqueEnhancingKnowledgeDiscovery2017|28327945|LOF|TBD||TRUE
|SF3B1|1|FALSE|NA||reddyGeneticFunctionalDrivers2017|28985567|NEO|23160465|cazzolaBiologicClinicalSignificance2013|FALSE
|SGK1|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119|GOF|33988691|gaoSGK1MutationsDLBCL2021|TRUE
|SIN3A|1|FALSE|NA||chapuyMolecularSubtypesDiffuse2018|29713087||||FALSE
|SMARCA4|1|FALSE|NA||reddyGeneticFunctionalDrivers2017|28985567|LOF|33144586|fernandoFunctionalCharacterizationSMARCA42020|FALSE
|SOCS1|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|15572583|melznerBiallelicMutationSOCS12005|TRUE
|SPEN|1|FALSE|NA||albuquerqueEnhancingKnowledgeDiscovery2017|28327945||||TRUE
|STAT3|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|GOF|23861822|huNovelMissenseM206K2013|TRUE
|STAT6|1|FALSE|NA||yildizActivatingSTAT6Mutations2015|25428220|GOF|35851155|mentzPARP14NovelTarget2022|TRUE
|TAF1|1|FALSE|NA||morinMutationalStructuralAnalysis2013|23699601||||FALSE
|TBL1XR1|1|FALSE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534|LOF|32619424|venturuttiTBL1XR1MutationsDrive2020|TRUE
|TET2|1|FALSE|NA||albuquerqueEnhancingKnowledgeDiscovery2017|28327945|LOF|23831920|asmarGenomewideProfilingIdentifies2013|TRUE
|TMEM30A|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|32094924|ennishiTMEM30ALossoffunctionMutations2020|TRUE
|TMSB4X|1|TRUE|NA||reddyGeneticFunctionalDrivers2017|28985567||||FALSE
|TNFAIP3|1|FALSE|NA||compagnoMutationsMultipleGenes2009|19412164|LOF|19412164|compagnoMutationsMultipleGenes2009|TRUE
|TNFRSF14|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|TBD||TRUE
|TOX|1|FALSE|NA||reddyGeneticFunctionalDrivers2017|28985567||||TRUE
|TP53|1|FALSE|NA||morinFrequentMutationHistonemodifying2011|21796119|LOF|12826609|katoUnderstandingFunctionstructureFunctionmutation2003|TRUE
|UBE2A|1|FALSE|NA||lohrDiscoveryPrioritizationSomatic2012|22343534||||TRUE
|WEE1|1|FALSE|NA||schmitzGeneticsPathogenesisDiffuse2018|29641966||||TRUE
|XPO1|1|FALSE|NA||mareschalWholeExomeSequencing2016|26608593|NEO|33007990|miloudiXPO1E571KMutationModifies2020|FALSE
|ZC3H12A|1|FALSE|NA||chapuyMolecularSubtypesDiffuse2018|29713087|LOF|19747262|skalniakRegulatoryFeedbackLoop2009|FALSE
|ZFP36L1|1|TRUE|NA||morinFrequentMutationHistonemodifying2011|21796119||||TRUE
|ZNF292|1|FALSE|NA||reddyGeneticFunctionalDrivers2017|28985567||||FALSE
|ZNF608|1|FALSE|NA||morinMutationalStructuralAnalysis2013|23699601||||FALSE
|LCOR|1|FALSE|NA||novakWholeexomeAnalysisReveals2015|26314988||||FALSE


### Burkitt lymphoma

The master curated list for BL including DLBCL genes that have been scrutinized within BL is can be found in `bl_genes.tsv` (or [here](bl_genes.tsv)). The earliest_support_BL column indicates the PubMed ID of the first study to nominate this gene as mutated in BL (or NA when not applicable). The frequency_BL_Thomas and frequency_BL_Panea columns report the percentage of patient samples with at least one non-silent mutation in this gene in the two studies. In the case of Panea et al, these numbers are based on the reanalysis of the exome data from this study by the Morin lab (as detailed in Dreval et al). The original frequencies based on the mutation calls from Panea et al are in the frequency_BL_Panea_original column.

Similar to dlbcl_genes.tsv, this file contains a Curated_BL_driver column, which will be TRUE for the subset of genes that have made it to our core list. The current core list at the time this document was created is below. _This may not match the actual file, depending on whether this document is kept up to date. Please refer to the the file rather than this list._ The links in the table below are not exhaustive but instead are there to show some examples of the mutations underlying this analysis. *Important:* The mutation frequency from our reanalysis of the data from Panea et al is reported under frequency_BL_Panea point whereas the lollipop plots are derived from the mutations originally reported by Panea. The different patterns in the mutations reported by that study and observed in our analysis has not been addressed by the authors. 

| Gene    | Curated_BL_driver | frequency_BL_Thomas | frequency_BL_Panea |
|---------|-------------------|---------------------|--------------------|
| MYC     | TRUE              | [60.2](../../literature/mutation_patterns/lollipop/by_study/as_reported/bl_panea/MYC.png)                 | [49.5](../../literature/mutation_patterns/lollipop/by_study/as_reported/bl_thomas/MYC.png)               |
| DDX3X   | TRUE              | 48.7                | 39.6               |
| ID3     | TRUE              | 47                  | 31.7               |
| TP53    | TRUE              | [41.9](../../literature/mutation_patterns/lollipop/by_study/as_reported/bl_thomas/TP53.png)                 | [43.6](../../literature/mutation_patterns/lollipop/by_study/as_reported/bl_panea/TP53.png)                |
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

### aSHM regions

We provide two (somewhat redundant) files that each represent the coordinates of genomic regions we have identified as being targets of aberrant somatic hypermutation (aSHM) from our analysis of DLBCL genomes. Use the bed file if you need that format, otherwise refer to the txt file for more complete details about the regions. _Contributors willing to help fill in regions you think need to be added to this list are encouraged to submit a GitHub issue and we'll review the data we have to determine if this is warranted._ 

The following three files are semi-redundant representations of the same information. This is the coordinates for regions of the genome commonly affected by aSHM in DLBCL. [This version](somatic_hypermutation_locations_with_DLBCL_frequencies.tsv) includes additional columns tracking the number and percentage of DLBCL genomes with mutations in each region based on our meta-analysis. 

`somatic_hypermutation_locations_GRCh37.txt`

`somatic_hypermutation_locations_with_DLBCL_frequencies.tsv`

`somatic_hypermutation_locations_GRCh37.bed`

## How to contribute

These lists may be incomplete and may contain errors. If you believe a gene is missing, is mis-attributed, etc please let us know. _Contributors willing to help fill in missing data (genes, citations, hot spots) are encouraged to submit a GitHub issue._ 

## How to cite

If you use the information in these lists please cite the following papers:

### DLBCL gene list

Minimum Information for Reporting a Genomics Experiment. Dreval K, Boutros PC, Morin RD.
Blood. 2022 Oct 11:blood.2022016095. doi: 10.1182/blood.2022016095. PMID: 36219881

### BL gene list

GENETIC SUBGROUPS INFORM ON PATHOBIOLOGY IN ADULT AND PEDIATRIC BURKITT LYMPHOMA. Thomas N, Dreval K et al. Blood. 2022 Oct 6:blood.2022016534. doi: 10.1182/blood.2022016534. PMID: 36201743

### SHM regions

Genome-wide discovery of somatic regulatory variants in diffuse large B-cell lymphoma. Arthur et al. Nat Commun. 2018 Oct 1;9(1):4001. doi: 10.1038/s41467-018-06354-3.

Super-enhancer hypermutation alters oncogene expression in B cell lymphoma. Bal E et al. Nature. 2022 Jul;607(7920):808-815. doi: 10.1038/s41586-022-04906-8. Epub 2022 Jul 6.
