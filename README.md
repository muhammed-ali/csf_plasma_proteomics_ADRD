# Large-scale CSF and plasma proteomics reveal dysregulation of immune system, synaptic impairment, and extracellular matrix related pathways in neurodegeneration

# Table of contents
* [Introduction](#introduction)
* [Content](#content)
* [Data](#data)
* [Requirements](#requirements)
* [License](#license)
* [Instructions](#instructions)

# Introduction
This repository contains the code for bioinformatics analyses described in the article "Large-scale CSF and plasma proteomics reveal dysregulation of immune system, synaptic impairment, and extracellular matrix related pathways in neurodegeneration".

This project investigated CSF and plasma proteomics data from the SomaScan assay to identify proteins associated with Alzheimer disease (AD), Parkinson's disease (PD), dementia with Lewy bodies (DLB), and Frontotemporal Dementia (FTD). Idnetified proteins were leveraged to characterize disease-specific and shared molecular signatures, and create disease-spcific prediction models. Biological pathway and cell type enrichment analyses were performed to understand underlying biology.

# Content
The code covers the following main analysis steps:

1. Data pre-processing: Proteomics data preparation and pricipal component analysis 
2. Differential expression analysis
3. Prediction model development using LASSO regression
4. Pathways and cell type enrichment analyses
   
# Data
Proteomics data analysed in this study is available at:
- ADNI: http://adni.loni.usc.edu/
- Knight-ADRC (CSF): https://dss.niagads.org/ (Accession: ng00130)
- Knight-ADRC (Plasma): https://live-knightadrc-washu.pantheonsite.io/professionals-clinicians/request-center-resources/
- FACE and Barcelona-1 cohorts: http://www.fundacioace.com/
- PPMI: https://www.ppmi-info.org/
- Stanford-ADRC: https://live-knightadrc-washu.pantheonsite.io/professionals-clinicians/request-center-resources/
- The harmonized GNPC dataset will be made publicly available following an embargo period at https://www.neuroproteome.org. 

# Requirements
The code was written in R (version 4.3.0) and relies on multiple R and Bioconductor packages, including:
- caret
- glmnet
- nlme
- pROC
- ROCR
- dplyr
- clusterProfiler
- ReactomePA
- ggplot2
- EnhancedVolcano

- Additional packages listed at the beginning of each R script

# License
The code is available under the MIT License.

# Instructions
The code was tested on R 4.3.0 on Linux operating systems, but should be compatible with later versions of R installed on current Linux, Mac, or Windows systems.

To run the code, the correct working directory containing the input data must be specified at the beginning of the R-scripts, otherwise the scripts can be run as-is.

The scripts should be run in the following order:

    data_preparation.R

    differential_expression_analysis.R

    prediction_models.R

    pathway_and_celltype_enrichment_analysis.R

