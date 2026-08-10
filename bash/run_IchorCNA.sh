#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH --time=08:00:00
#SBATCH --output=logs/ichor%j.out

sample=""
home_dir=""
R_script_paths=""
reference_path=""
ichor_CNA_path=""

usage() {
    echo "Usage: $0 Required arguments: -s SAMPLE_NAME -o OUTPUT_DIRECTORY -r REFERENCE_PATH -R R_PATH -I ICHORCNA_PATH"
    echo
    echo "Options:"
    echo "  -s SAMPLE   Sample name matching a BAM file alignment/<SAMPLE>_rmdup.bam"
    echo "  -o DIRECTORY   Output directory"
    echo "  -r PATH   Path to reference data (https://github.com/CTrierMaansson/cfChIP_pembro/reference/)"
    echo "  -R PATH   Path to R scripts (https://github.com/CTrierMaansson/cfChIP_pembro/R/scripts/)"
    echo "  -I PATH   Path to IchorCNA directory (https://github.com/broadinstitute/ichorCNA)"
}

while getopts ":s:o:r:R:I:h" opt; do
    case "$opt" in
        s)
            sample="$OPTARG"
            ;;
        o)
            home_dir="$OPTARG"
            ;;
        r)
            reference_path="$OPTARG"
            ;;
        R)
            R_script_paths="$OPTARG"
            ;;
        I)
            ichor_CNA_path="$OPTARG"
            ;;

        h)
            usage
            exit 0
            ;;
        :)
            echo "Error: option -$OPTARG requires a value" >&2
            exit 1
            ;;
        \?)
            echo "Error: unknown option -$OPTARG" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$sample" ||
      -z "$home_dir" ||
      -z "$reference_path" ||
      -z "$R_script_paths" ||
      -z "$ichor_CNA_path" ]]; then
    echo "Usage: $0 Required arguments: -s SAMPLE_NAME -o OUTPUT_DIRECTORY -r REFERENCE_PATH -R R_PATH -I ICHORCNA_PATH"
    exit 1
fi

bam_file="${sample}_rmdup.bam"

require_dir() {
    local path="$1"
    local option="$2"

    if [[ ! -d "$path" ]]; then
        echo "Error: $option must point to an existing directory: $path" >&2
        exit 1
    fi
}

echo "Testing whether input files can be found"
echo "cfChIP BAM file: $home_dir/alignment/$bam_file"
echo "Reference chromosome sizes: ${reference_path}/hg38_size.bed"
echo "Normal panel: ${reference_path}/NormalPanel.rds"
echo "Background regions R script: ${R_script_paths}/individual_on_targets.R"
echo "IchorCNA directory: ${ichor_CNA_path}"

require_dir "$home_dir" "-o"
require_dir "$R_script_paths" "-R"
require_dir "$ichor_CNA_path" "-I"
require_dir "$reference_path" "-r"

if [[ ! -f "$home_dir/alignment/$bam_file" ]]; then
    echo "Error: cannot find $bam_file in $home_dir/alignment/" >&2
    exit 1
fi

#The R scripts are available at:
#https://github.com/CTrierMaansson/cfChIP_pembro/R/scripts/

if [[ ! -f "$R_script_paths/individual_on_targets.R" ]]; then
    echo "Error: individual_on_targets.R cannot be found in: $R_script_paths" >&2
    exit 1
fi

#The hg38 chromosome sizes used in this project are available at:
#https://github.com/CTrierMaansson/cfChIP_pembro/reference/

if [[ ! -f "$reference_path/hg38_size.bed" ]]; then
    echo "Error: hg38_size.bed cannot be found in: $reference_path" >&2
    exit 1
fi

#We created a panel of normals to use for normalization
#for cancer cases. This panel can be found at:
#https://github.com/CTrierMaansson/cfChIP_pembro/reference/

if [[ ! -f "$reference_path/NormalPanel.rds" ]]; then
    echo "Error: NormalPanel.rds cannot be found in: $reference_path" >&2
    exit 1
fi

echo "Test successful"

echo "This script will prepare files for IchorCNA analysis of a cfChIP sample: $sample"

echo "Finding peak regions and defining background regions"

eval "$(conda shell.bash hook)"
conda activate cfchip_ichor

mkdir -p $home_dir/IchorCNA/background_regions/
mkdir -p $home_dir/IchorCNA/wig/
mkdir -p $home_dir/IchorCNA/results/


Rscript \
    --vanilla \
    $R_script_paths/individual_on_targets.R \
    $home_dir/alignment/$bam_file \
    $reference_path/hg38_size.bed

echo "Extracting reads in background regions"

background_bed="${sample}_background_regions.bed"
output_bam="${sample}_rmdup_background.bam"

echo "BAM file with background reads: $output_bam"

samtools view \
  -b \
  -L $home_dir/IchorCNA/background_regions/$background_bed \
  $home_dir/alignment/$bam_file \ 
  > $home_dir/alignment/$output_bam
  
echo "indexing BAM file"

samtools index \
    $home_dir/alignment/$output_bam

output_wig="${sample}_background.wig"

echo "Preparing 50kb bin background WIG file: $output_wig"

readCounter \
  --window 50000 \
  --quality 1 \
  --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
   $home_dir/alignment/$output_bam \
   > $home_dir/IchorCNA/wig/$output_wig

echo "Running IchorCNA"

Rscript $ichor_CNA_path/scripts/runIchorCNA.R \
  --id $sample \
  --WIG $home_dir/IchorCNA/wig/$output_wig \
  --normalPanel $reference_path/NormalPanel.rds \
  --gcWig $ichor_CNA_path/inst/extdata/gc_hg38_50kb.wig \
  --mapWig $ichor_CNA_path/inst/extdata/map_hg38_50kb.wig \
  --centromere $ichor_CNA_path/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
  --includeHOMD FALSE \
  --estimateScPrevalence FALSE \
  --scStates "c()" \
  --ploidy "c(2)" \
  --normal "c(0.90,0.95, 0.99, 0.995, 0.999)" \
  --maxCN 3 \
  --genomeBuild "hg38" \
  --genomeStyle "UCSC" \
  --outDir $home_dir/IchorCNA/results
  
echo "IchorCNA analysis DONE"