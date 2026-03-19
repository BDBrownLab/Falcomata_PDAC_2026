# Visium Spatial Transcriptomics Analysis

This repository contains the workflows used to pre-process and analysis public 10X Visium Spatial Transcriptomics data for Falcomata et al, 2026. The following datasets were downloaded for this analysis:

| Dataset        | Publication | 
| ------------- |:-------------:| 
| HTAN and Zhou et al, 2022 | PMID: 35902743 |  
| Pei et al, 2025 | PMID: 40269162 |  
| Chen et al, 2025 | PMID: 40680743|  

Examples of the pre-processing, integration, and analysis frameworks can be found in the following files. The majority of the analyses were conducted using pre-existing packages, which are cited when used.

| Analysis        | Script | Notes |
| ------------- |:-------------:| :-------------:| 
| Pre-processing | 01-visium_preprocessing.ipynb | Largely follows Squidpy. Keep raw data for cell2location.
| Dataset integration and deconvolution | 02-visium_integration_decon.ipynb | Largely adapted from the SCVI workflow and cell2location workflow
| Analysis | 03-visium_paper_analysis.ipynb|  

