#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH -c 1
#SBATCH --output=logs/peak_call%j.out
#SBATCH --time=08:00:00
#SBATCH --array=0-59


files=("Donor1" "Donor2" "Donor3" "Donor4" "Donor5" "Donor6" "pt03_BL" "pt03_w6" "pt04_BL" "pt04_w6" "pt06_BL" "pt06_w9" "pt08_BL" "pt08_w9" "pt09_BL" "pt09_w6" "pt11_BL" "pt11_w6" "pt12_BL" "pt12_w6" "pt16_BL" "pt16_w6" "pt19_BL" "pt19_w6" "pt21_BL" "pt21_w3" "pt22_BL" "pt22_w6" "pt23_BL" "pt23_w9" "pt24_BL" "pt24_w6" "pt25_BL" "pt25_w6" "pt26_BL" "pt26_w6" "pt28_BL" "pt28_w9" "pt29_BL" "pt29_w6" "pt31_BL" "pt31_w9" "pt33_BL" "pt33_w6" "pt36_BL" "pt36_w6" "pt37_BL" "pt37_w6" "pt38_BL" "pt38_w6" "pt39_BL" "pt39_w6" "pt41_BL" "pt41_w3" "pt45_BL" "pt45_w6" "pt46_BL" "pt46_w6" "pt49_BL" "pt49_w6")

home_dir="/define/home/directory"

sample="${files[$SLURM_ARRAY_TASK_ID]}"

echo "# This is the peak calling and annotation of the following sample: ${sample}"

eval "$(conda shell.bash hook)"
conda activate cfChIP_pembro

mkdir -p $home_dir/macs
mkdir -p $home_dir/macs/$sample
mkdir -p $home_dir/macs/annotated

cfChIP_file="${sample}_rmdup.bam"

echo "This is the peak calling of ${cfChIP_file}"

macs2 callpeak \
    -t $home_dir/alignment/$cfChIP_file \
    -g hs \
    --keep-dup all \
    --outdir $home_dir/macs/$sample \
    -n $sample \
    -B \
    -q 0.1 \
    --broad 

peak_file="${sample}_peaks.broadPeak"

echo "Annotating peaks"

#The R scripts are available at:
#https://github.com/CTrierMaansson/cfChIP_pembro/R/scripts/
R_script_paths="path/to/R_scripts"

Rscript 
    --vanilla 
    $R_script_paths/annotate_peaks.R \
    $home_dir/macs/$sample/$peak_file
