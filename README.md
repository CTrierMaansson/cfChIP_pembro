# cfChIP_pembro

This repository contains the analysis scripts and reference files used in the study “Plasma histone modification profiling can infer transcriptional programs during immune checkpoint inhibitor therapy in non-small cell lung cancer.” The scripts implement the cfChIP-seq preprocessing, peak calling, copy-number analysis, and downstream statistical analyses described in the manuscript.

## Table of contents
- [Repository content](#Repository-content)
- [Environments](#Environments)
- [Code](#Code)
- [Result directory structure](#Result-directory-structure)

## Repository structure

This repository contains the scripts and reference data used to generate the
results for the manuscripts. 

 - [bash/](bash/) - Scripts to generate data used for downstream analyses
 - [reference/](reference/) - Reference data generated or downloaded for the analyses 
 - [R/](R/) - R scripts and commonly used functions for analysis or to generate plots
 - [meta.rds](meta.rds) - metadata for cfChIP files and ctDNA information for samples

## Requirements and environments
The required software dependencies are specified in the corresponding 
conda environment files:

1. cfchip_pembro.yml - Main environment for raw data processing
2. cfchip_ichor.yml - Environment to run IchorCNA on cfChIP-seq background reads

These environments can be installed using 

```bash
conda env create \
    -n cfchip_pembro \
    -f /path/to/cfChIP_pembro/cfchip_pembro.yml

conda env create \
    -n cfchip_ichor \
    -f /path/to/cfChIP_pembro/cfchip_ichor.yml
```

The required R dependencies for downstream analysis
are listed within each file in [R/](R/)

## Analysis workflow

### Preprocessing

Three bash scripts are used to generate the files used for downstream analyses.

The first script, preprocessing.sh performs QC, alignment, filtering,
and creates files needed for downstream analyses.


```bash
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

### Peak calling

Peak calling is then performed using:

```bash
sbatch /cfChIP_pembro/bash/peak_call.sh \
  -s <sample_name> \
  -o /output/path/ \
  -R /path/to/cfChIP_pembro/R/scripts/ \
```
### Copy-number analysis

Finally, copy-number alterations are analyzed using IchorCNA on background reads

```bash
sbatch /cfChIP_pembro/bash/run_IchorCNA.sh \
  -s <sample_name> \
  -o /output/path/ \
  -r /path/to/cfChIP_pembro/reference/ \
  -R /path/to/cfChIP_pembro/R/scripts/ \
  -I /path/to/IchorCNA/
```

### Downstream analysis

[R/scripts/](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/R/scripts/) contains helper R scripts called by the bash pipelines.

In addition, [R/](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/R/) 
contains three files containing the functions used for downstream data analysis:

 - [admin_functions.R](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/R/admin_functions.R) - Definitions of variables and data formats
 - [analysis_functions.R](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/R/analysis_functions.R) - Downstream analysis of data generated with preprocessing.sh
 - [plot_functions.R](https://github.com/CTrierMaansson/cfChIP_pembro/blob/main/R/plot_functions.R) - Creation of plots based on results from analysis_functions.R
 
## Output structure

Example of output structure for a single sample (D1) where the output 
path is /cfChIP_pembro/res/

```bash
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
    ├── fragment_lengths (fragmentlength distributions)
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

## Data availability
Raw sequencing data are not included in this repository. We are in process of 
publishing the raw genome coverage profiles to [EGA](https://ega-archive.org/).
The repository contains the analysis code and reference files required to 
reproduce the processing and downstream analyses.

## Citation

If you use this repository or code in your research, please cite:

Maansson CT, et al. Plasma histone modification profiling can infer transcriptional programs during immune checkpoint inhibitor therapy in non-small cell lung cancer. 2026.


