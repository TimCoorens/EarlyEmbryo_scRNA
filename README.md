This repository contains the code used for the analysis of the single-cell RNA sequencing data of early embryos, used in the paper 'Single Cell Analyses of Human Embryos Defines the Putative Anterior Hypoblast Signalling Centre'.
For any queries, please email tc16@sanger.ac.uk.

The analysis of scRNAseq data in these scripts relies on Seurat (https://satijalab.org/seurat/). 

The script 'scRNA_embryo_data_process.R' contains all code used for integration and batch correction of cellranger output on our own data. It also has the code to generate Fig. 1a-f, extended data fig. 1

The script 'scRNA_embryo_expression_patterns.R' contains all code used to assess expression levels of certain genes (i.e. the FGF family), as well a closer inspection into subclustering of hypo- and epiblast. It  has the code to generate Fig. 2a-b, Fig. 3a-b,i, and extended data fig. 6.

The script 'scRNA_embryo_logistic_regression.R' contains all code used for comparison between our own data and other datasets, as well as the integration with data from Stirparo et al. It  has the code to generate Fig 1g, extended data fig. 3, 5.


