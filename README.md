# Epigenomics of Homologous Recombination Deficiency (HRD) in Triple-Negative Breast Cancer (TNBC)

Y. David Chen, Ph.D. 

Manuscript title: Chen et al. Extensive epigenomic dysregulation is a hallmark of homologous recombination deficiency in triple-negative breast cancer

## Repository Structure 

* `dataproc/` Essential bioinformatic data preprocessing procedures e.g. MethylationEPIC processing
* `main/` Core statistical analyses and associated data visualization, e.g. differential methylation analysis
* `suppl/` Supplemental analyses, including machine learning (ML) modeling
* `utils.R` and `plotThemes.R` contain reusable, global helpers and constants for importing in any R scripts

Scripts used for data download/acquisition, data cleansing, exploratory data analysis (EDA) not included in the final paper, and constants/PATHs/CONFIGs storage may be excluded/masked.

## List of Analyses

### Bioinformatics Processing (`dataproc/`)

* MethylationEPIC QC, preprocessing, and data prep
* Preprocessing and extraction of segmented copy number from methylation arrays

### Core Analyses (`main/`)

* Study population summarization (Table One)
* Fraction of Genome Altered
* Methylation Dysregulation Index
* Promoter hypermethylation of _BRCA1_ and _RAD51C_
* Global methylation clusters with RPMM
* Methylation at transcriptional regulatory regions
* Age-adjusted EWAS
* Genomic context enrichment
* Methylation heatmap visualization at various levels

### Supplemental Analyses (`suppl/`)

* Supervised prediction of subject race with SNP probes
* Supervised learning of continuous HRD probability scores and biological discovery on TCGA

## Disclaimers

Materials from this repository may be subject to copyright. Use the materials responsibly and at your own risk. 

If you referenced the datasets, findings, or methods including code from this work, please include the appropriate citation(s).


