# cfChIP_pembro

This repository contains the fundamental scripts associated with the manuscript entitled
"Plasma nucleosome profiling to infer transcriptional dynamics during immunotherapy in lung cancer"
By Christoffer Trier Maansson et al. (2026)

## Table of contents
- [Environments](#Environments)
- [Content](#Content)
- [Running code](#Running code)

## Environments
The data is analyzed using two conda environments
1. cfchip_pembro - Main environment for raw data processing
2. cfchip_ichor - Environment to run IchorCNA on cfChIP-seq background reads

These environments can be installed using 

```{bash}
conda env create \
    -n cfchip_pembro \
    -f /cfChIP_pembro/cfchip_pembro.yml

conda env create \
    -n cfchip_ichor \
    -f /cfChIP_pembro/cfchip_ichor.yml
```

## Content

This repository contains the scripts and reference data used to generate the
results for the manuscripts. 
 - /cfChIP_pembro/bash/ - Scripts to generate data used for downstream analyses
 - /cfChIP_pembro/reference/ - Reference data generated for the analyses or accessed publically
 - /cfChIP_pembro/R/ - R scripts and commonly used functions used to generate plots
 
### Reference data
Here is some info about the files in reference because they require more info

## Running code

Three bash scripts will create the files used for downstream analyses.

Something more about this when arguments have been implemented

```{bash}
sbatch /cfChIP_pembro/bash/preprocessing.sh


```

### result directory structure

Here I want the tree of the result directory



