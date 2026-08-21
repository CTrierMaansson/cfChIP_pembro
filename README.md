# cfChIP_pembro

This repository contains the fundamental scripts associated with the manuscript entitled
"Plasma histone modification profiling can infer transcriptional programs during immune checkpoint inhibitor therapy in non-small cell lung cancer"
By Christoffer Trier Maansson et al. (2026)

## Table of contents
- [Environments](#Environments)
- [Repository content](#Repository-content)
- [Code](#Code)
- [Result directory structure](#Result-directory-structure)

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

## Repository content

This repository contains the scripts and reference data used to generate the
results for the manuscripts. 
 - [bash/](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/bash/) - Scripts to generate data used for downstream analyses
 - [reference/](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/reference/) - Reference data generated or downloaded for the analyses 
 - [R/](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/R/) - R scripts and commonly used functions for analysis or to generate plots
 - [meta.rds](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/meta.rds) - metadata for cfChIP files and ctDNA information for samples

## Code

### bash scripts

Three bash scripts will create the files used for downstream analyses.

The first script, preprocessing.sh performs QC, alignment, filtering,
and creates files needed for downstream analyses.


```{bash}
sbatch /cfChIP_pembro/bash/preprocessing.sh \
  -s <sample_name> \
  -o /output/path/ \
  -f /path/to/<sample_name>_R1.fastq \
  -F /path/to/<sample_name>_R2.fastq \
  -r /path/to/cfChIP_pembro/reference/ \
  -R /path/to/cfChIP_pembro/R/scripts/ \
  -g /path/to/indexed_reference_genome/ \
  -b /path/to/fgbio.jar \
  -m /path/to/hmftools.jar \
  -p hg38

```

Following this, the peak calling is performed using:

```{bash}
sbatch /cfChIP_pembro/bash/peak_call.sh \
  -s <sample_name> \
  -o /output/path/ \
  -R /path/to/cfChIP_pembro/R/scripts/ \
```

Lastly, Copy number alterations are analyzed using IchorCNA on background reads

```{bash}
sbatch /cfChIP_pembro/bash/run_IchorCNA.sh \
  -s <sample_name> \
  -o /output/path/ \
  -r /path/to/cfChIP_pembro/reference/ \
  -R /path/to/cfChIP_pembro/R/scripts/ \
  -I /path/to/IchorCNA/
```

### R

[R/scripts/](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/R/scripts/) contain the R scripts used in the bash scripts above.

In addition, [R/](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/R/) 
contains three files containing the functions used for downstream data analysis:

 - [admin_functions.R](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/R/admin_functions.R) - Definitions of variables and data formats
 - [analysis_functions.R](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/R/analysis_functions.R) - Downstream analysis of data generated with preprocessing.sh
 - [plot_functions.R](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/R/plot_functions.R) - Creation of plots based results from analysis_functions.R

## Result directory structure

Example of output structure for a single sample (D1) where the output 
path is /cfChIP_pembro/res/

```{bash}
/cfChIP_pembro/
├── bash
├── logs
├── R
│   └── scripts
├── reference
└── res
    ├── alignment (.sam and .bam files)
    ├── bed (fragment .bed files)
    ├── bedgraph (genome coverage profiles)
    ├── bigwig (genome coverage profiles)
    ├── coverage (Coverage relative to TSS)
    ├── enrichment_ratios (sample enrichment ratios)
    ├── fastqc (output of fastqc)
    ├── fragment_lengths (fragmentlength distrivutions)
    ├── IchorCNA (CNA results from IchorCNA)
    │   ├── background_regions
    │   ├── results
    │   │   └── D1
    │   └── wig
    ├── macs (peak calling results)
    │   ├── annotated
    │   └── D1
    ├── reduced_fragments (centered fragment .bed files)
    ├── trimmed_fastqs (trimmed .fastq files)
    │   └── fastqc
    └── umi_stats (stats from UMI deduplication)

```



