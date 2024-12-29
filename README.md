# Enhance ctDNA quantification by ctCandi score and ichorCNA tumor fraction

## List of datasets

- `LOD` dataset: diluted cfDNA samples with known tumor abundance (Folder name: `reads_from_450_regions_LOD_samples`). 

- `REPORT 4 metadata` dataset: dataset used in the project ECD, REPORT 4 evaluation (Folder name: `reads_from_450_regions_with_readname_TMDfull`). 

- `in-silico spike-in` dataset: tumor reads were sampled from tumor tissue samples and spike-in healthy control cfDNA samples at multiple fractions (Folder name: `reads_from_450_regions_spikein_v2`). 

- `validation` dataset: a set of samples from MRD and ECD projects for final validation (Folder name: `reads_from_450_regions_validations_Vi-Truong`). 

## Analysis pipeline

- `01-06` analysis scripts: pre-process and find differential methylated loci/regions in TCGA data and GS-Target data. 

- `07_count_reads_in_regions.py`: 