# hiPSC-derived_AxD_models
This repository contains code used to process and analyze two scRNA-seq datasets prepared from co-culture and organoid models of Alexander disease. Samples were processed using the 10X Genomics platform. The data have been deposited in Gene Expression Omnibus under accession number GSE261158. The results have been published as **paper**.

The repository includes the following folders and files:

### AxD_cocultures
_samples.csv_
- information about samples which was used in preprocessing

_2017_07_preprocessing.R_
- fastq_screen
- STAR alignment

_20211222_scRNA_preprocessing.Rmd_
- EmptyDrops - filtering of empty droplets
- quality control
- integration
- initial clustering and annotation
- identification of low-quality cells
- DoubletFinder

_20211230_scRNA_analysis.Rmd_
- more accurate annotation
- subsetting to cell populations
- velocyto, monocle

_20220221_coculture_AxD_code.Rmd_
- identification of contamination with SoupX and repeated analysis of corrected data including quality control
- clustering, subclustering, and annotation

_20220307_scRNA_post_SoupX.Rmd_
- Cerebro
- additional quality control and finalized annotation of clusters

_20220225_DEA.Rmd_
- identification of differentially expressed genes 
- Gene Ontology enrichment analysis

_20230216_CellChat.Rmd_
- CellChat analysis of cell-cell interaction patterns

_20221006_coculture_figures.Rmd_
- code used to generate all figures that are in the manuscript

### AxD_organoids
_samples.csv_
- information about samples which was used in preprocessing

_2017_07_preprocessing.R_
- fastq_screen
- STAR alignment

_scRNA_analysis.Rmd_
- Seurat analysis
- integration
- quality control
- clustering
- annotation
- comparison of RNAseq and proteomics data
- code used to generate figures in the manuscript

_DEA.R_
- differential expression analysis
- Gene Ontology enrichment analysis

_RTqPCR_plots.R_
- code used to generate RT-qPCR plots

_Gruffi_stress.Rmd_
- Gruffi tool used to identify stressed cells in organoids

_CellChat.Rmd_
- CellChat analysis of cell-cell interaction patterns

_overlap_with_coculture.R_
- code used to compare co-culture and organoid datasets





