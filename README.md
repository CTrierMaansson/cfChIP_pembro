# cfChIP_pembro

This repository contains the analysis scripts and reference files used in the study “Plasma histone modification profiling can infer transcriptional programs during immune checkpoint inhibitor therapy in non-small cell lung cancer.” The scripts implement the cfChIP-seq preprocessing, peak calling, copy-number analysis, and downstream statistical analyses described in the manuscript.

## Table of contents
- [Repository structure](#Repository-structure)
- [Requirements and environments](#Requirements-and-environments)
- [Analysis workflow](#Analysis-workflow)
- [Output structure](#Output-structure)
- [Data availability](#Data-availability)
- [Citation](#Citation)

## Repository structure

This repository contains the scripts and reference data used to generate the
results for the manuscript. 

 - [bash/](bash/) - Scripts to generate data used for downstream analyses
 - [reference/](reference/) - Reference data generated or downloaded for the analyses 
 - [R/](R/) - R scripts and commonly used functions for analysis or to generate plots
 - [`meta.rds`](meta.rds) - metadata for cfChIP-seq files and ctDNA information for samples

## Requirements and environments
The required software dependencies are specified in the corresponding 
conda environment files:

1. `cfchip_pembro.yml` - Main environment for raw data processing
2. `cfchip_ichor.yml` - Environment to run IchorCNA on cfChIP-seq background reads

These environments can be installed using 

```bash
conda env create \
    -n cfchip_pembro \
    -f /path/to/cfChIP_pembro/cfchip_pembro.yml

conda env create \
    -n cfchip_ichor \
    -f /path/to/cfChIP_pembro/cfchip_ichor.yml
```

In addition to the conda environments, the preprocessing workflow requires 
an indexed [hg38 reference genome](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/), 
[fgbio](https://github.com/fulcrumgenomics/fgbio/releases), 
and [HMFtools](https://github.com/hartwigmedical/hmftools/releases/tag/redux-v2.0).
[IchorCNA](https://github.com/broadinstitute/ichorCNA) is required for copy-number analysis. 


Required R packages for downstream analyses are specified in the
corresponding R scripts in [R/](R/)

## Analysis workflow

The analysis workflow consists of three primary processing steps: 
(1) preprocessing of raw cfChIP-seq FASTQ files, including QC, alignment and filtering 
(2) peak calling 
(3) copy-number analysis using IchorCNA on background reads. 

Subsequently, downstream statistical analyses and visualization 
are performed using R.

### Preprocessing

The [`preprocessing.sh`](bash/preprocessing.sh) script performs QC, alignment, filtering,
and other preprocessing steps required for downstream analyses:

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

Peak calling is then performed using [`peak_call.sh`](bash/peak_call.sh):

```bash
sbatch /cfChIP_pembro/bash/peak_call.sh \
  -s <sample_name> \
  -o /output/path/ \
  -R /path/to/cfChIP_pembro/R/scripts/ \
```
### Copy-number analysis

Finally, copy-number alterations are analyzed using IchorCNA on background
reads with [`run_IchorCNA.sh`](bash/run_IchorCNA.sh):

```bash
sbatch /cfChIP_pembro/bash/run_IchorCNA.sh \
  -s <sample_name> \
  -o /output/path/ \
  -r /path/to/cfChIP_pembro/reference/ \
  -R /path/to/cfChIP_pembro/R/scripts/ \
  -I /path/to/IchorCNA/
```

### Downstream analysis

[R/scripts/](R/scripts/) contains helper R scripts called by the bash pipelines.

In addition, [R/](R/) 
contains three files containing the functions used for downstream data analysis:

 - [`admin_functions.R`](R/admin_functions.R) - Definitions of variables and data formats
 - [`analysis_functions.R`](R/analysis_functions.R) - Downstream analysis of data generated with preprocessing.sh
 - [`plot_functions.R`](R/plot_functions.R) - Creation of plots based on results from analysis_functions.R
 
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
    ├── alignment   # SAM/BAM alignment files
    ├── bed     # Fragment BED files
    ├── bedgraph    # Genome coverage profiles
    ├── bigwig  # Genome coverage profiles
    ├── coverage    # Coverage relative to TSS
    ├── enrichment_ratios   #Sample enrichment ratios
    ├── fastqc  # Output of fastqc
    ├── fragment_lengths    # Fragment length distributions
    ├── IchorCNA    # CNA results from IchorCNA
    │   ├── background_regions
    │   ├── results
    │   │   └── D1
    │   └── wig
    ├── macs    # Peak calling results
    │   ├── annotated
    │   └── D1
    ├── reduced_fragments   # Centered fragment .bed files
    ├── trimmed_fastqs  #Trimmed .fastq files
    │   └── fastqc
    └── umi_stats   #Stats from UMI deduplication

```

## Data availability
Raw sequencing data are not included in this repository. Sequencing data are available through the  [European Genome-phenome Archive](https://ega-archive.org/) (EGA; accession: XXXX).
The repository contains the analysis code and reference files required to 
reproduce the processing and downstream analyses.

## Citation

If you use this repository or code in your research, please cite:

Maansson CT, et al. Plasma histone modification profiling can infer transcriptional programs during immune checkpoint inhibitor therapy in non-small cell lung cancer. 2026.(DOI: XXXX)


