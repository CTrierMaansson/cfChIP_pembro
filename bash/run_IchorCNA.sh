#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 1
#SBATCH --time=04:00:00
#SBATCH --array=0-59
#SBATCH --output=logs/ichor%j.out
#SBATCH --job-name="ichor"

files=("Donor1" "Donor2" "Donor3" "Donor4" "Donor5" "Donor6" "pt03_BL" "pt03_w6" "pt04_BL" "pt04_w6" "pt06_BL" "pt06_w9" "pt08_BL" "pt08_w9" "pt09_BL" "pt09_w6" "pt11_BL" "pt11_w6" "pt12_BL" "pt12_w6" "pt16_BL" "pt16_BL_rep2" "pt16_w6" "pt19_BL" "pt19_w6" "pt21_BL" "pt21_w3" "pt22_BL" "pt22_w6" "pt23_BL" "pt23_w9" "pt24_BL" "pt24_w6" "pt25_BL" "pt25_w6" "pt26_BL" "pt26_w6" "pt28_BL" "pt28_w9" "pt29_BL" "pt29_w6" "pt31_BL" "pt31_w9" "pt33_BL" "pt33_BL_rep2" "pt33_w6" "pt33_w6_rep2" "pt36_BL" "pt36_w6" "pt37_BL" "pt37_w6" "pt38_BL" "pt38_w6" "pt39_BL" "pt39_BL_rep2" "pt39_w6" "pt41_BL" "pt41_w3" "pt45_BL" "pt45_w6" "pt46_BL" "pt46_w6" "pt49_BL" "pt49_w6")

eval "$(conda shell.bash hook)"
conda activate cfchip_ichor

sample="${files[$SLURM_ARRAY_TASK_ID]}"

home_dir="/define/home/directory"

echo "This script will prepare files for IchorCNA analysis of a cfChIP sample: $sample"

echo "Finding peak regions and defining background regions"

mkdir -p $home_dir/IchorCNA/background_regions/
mkdir -p $home_dir/IchorCNA/wig/
mkdir -p $home_dir/IchorCNA/results/

#The R scripts are available at:
#https://github.com/CTrierMaansson/cfChIP_pembro/R/scripts/
R_script_paths="path/to/R_scripts"
bam_file="${sample}_rmdup.bam"

#The hg38 chromosome sizes used in this project are available at:
#https://github.com/CTrierMaansson/cfChIP_pembro/reference/
reference_path="path/to/reference/"

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

#IchorCNA was cloned from:
#https://github.com/broadinstitute/ichorCNA
ichor_CNA_path="/path/to/ichorCNA/"

#We created a panel of normals to use for normalization
#for cancer cases. This panel can be found at:
#https://github.com/CTrierMaansson/cfChIP_pembro/reference/

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