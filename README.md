# Enhanced scRNA-seq Workflow used in Ma et al. 2020
# Single-cell RNA-seq analysis of human coronary arteries using an enhanced workflow reveals SMC transitions and candidate drug targets
Comprehensive single-cell RNA seq analysis workflow: automatic cell labeling to drug targeting 
Wei Feng Ma, Chani J. Hodonsky, Adam W. Turner, Doris Wong, Yipei Song, Nelson Barrientos, Clint L. Miller

University of Virginia
Center for Public Health Genomics

## Overview
This workflow is divided into TWO parts:
1) data import into Seurat, automated cell labeling in singleR, Trajectory analysis in Monocle3
2) ligand-receptor profiling in scTalk, drug target analysis in DGIdb
![](images/scRNA_workflow.png)

## Data Import
If you are planning to use this pipeline, we assume that the data is coming from cellranger or from GEO as a data matrix. Some data sources come as .h5 or as a seurat object (.rds) that you can directly import, and you can skip this step. 

## Support
Raise an issue in our github and we will get back to you ASAP!
