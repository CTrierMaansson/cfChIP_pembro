#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 32G
#SBATCH -c 1
#SBATCH --output=logs/preprocess%j.out
#SBATCH --time=36:00:00
#SBATCH --array=0-59


files=("Donor1" "Donor2" "Donor3" "Donor4" "Donor5" "Donor6" "pt03_BL" "pt03_w6" "pt04_BL" "pt04_w6" "pt06_BL" "pt06_w9" "pt08_BL" "pt08_w9" "pt09_BL" "pt09_w6" "pt11_BL" "pt11_w6" "pt12_BL" "pt12_w6" "pt16_BL" "pt16_w6" "pt19_BL" "pt19_w6" "pt21_BL" "pt21_w3" "pt22_BL" "pt22_w6" "pt23_BL" "pt23_w9" "pt24_BL" "pt24_w6" "pt25_BL" "pt25_w6" "pt26_BL" "pt26_w6" "pt28_BL" "pt28_w9" "pt29_BL" "pt29_w6" "pt31_BL" "pt31_w9" "pt33_BL" "pt33_w6" "pt36_BL" "pt36_w6" "pt37_BL" "pt37_w6" "pt38_BL" "pt38_w6" "pt39_BL" "pt39_w6" "pt41_BL" "pt41_w3" "pt45_BL" "pt45_w6" "pt46_BL" "pt46_w6" "pt49_BL" "pt49_w6")

home_dir="/define/home/directory"

sample="${files[$SLURM_ARRAY_TASK_ID]}"

echo "# This is the preprocssing of the following sample: ${sample}"

eval "$(conda shell.bash hook)"
conda activate cfChIP_pembro

mkdir -p $home_dir/alignment
mkdir -p $home_dir/fastqc
mkdir -p $home_dir/trimmed_fastqs/fastqc
mkdir -p $home_dir/umi_stats
mkdir -p $home_dir/bed
mkdir -p $home_dir/reduced_fragments
mkdir -p $home_dir/enrichment_ratios
mkdir -p $home_dir/bigwig
mkdir -p $home_dir/bedgraph
mkdir -p $home_dir/fragment_lengths

input_R1_file="${sample}_R1.fastq"
input_R2_file="${sample}_R2.fastq"

echo "## Running pre trimming fastQC"
fastq_path="path/to/fastq"

fastqc
    -o $home_dir/fastqc/ \
    -f fastq \
    $fastq_path/$input_R1_file \
    $fastq_path/$input_R2_file

echo "## Running fastq trimming"

trim_galore \
    --fastqc \
    --fastqc_args "--outdir ${home_dir}/trimmed_fastqs/fastqc/" \
    --paired \
    -o trimmed_fastqs \
    $fastq_path/$input_R1_file \
    $fastq_path/$input_R2_file

trimmed_R1_file="${sample}_R1_val_1.fq"
trimmed_R2_file="${sample}_R2_val_2.fq"
sam_file="${sample}.sam"
reference_path="path/to/reference"

bowtie2 \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -x $reference_path \
    -p 4 \
    -1 $home_dir/trimmed_fastqs/$trimmed_R1_file \
    -2 $home_dir/trimmed_fastqs/$trimmed_R2_file \
    -S $home_dir/alignment/$sam_file 

bam_file="${sample}.bam"

echo "## Creating BAM file"

white_list="path/to/whitelist.bed"
#The hg38 whitelist used for this manuscript is available at
#https://github.com/CTrierMaansson/cfChIP_pembro/reference/

samtools view \
    -b \
    -L $white_list \
    -o $home_dir/alignment/$bam_file \
    $home_dir/alignment/$sam_file

echo "## BAM file done"
echo "## Sorting BAM file"

sorted_file="${sample}_sorted.bam"

samtools sort \
    -o $home_dir/alignment/$sorted_file \
    $home_dir/alignment/$bam_file

echo "## Copying UMIs"

annotate_file="${sample}_sorted_annotated.bam"
fgbio_path="path/to/fgbio.jar"
#We downloaded fgbio from https://github.com/fulcrumgenomics/fgbio/releases

java \ 
    -Xmx32g \
    -XX:+AggressiveHeap\
    -jar $fgbio_path \
    CopyUmiFromReadName \
    -i $home_dir/alignment/$sorted_file \
    -o $home_dir/alignment/$annotate_file 

echo "## Removing secondary alignments"

filtered_file="${sample}_sorted_annotated_filtered.bam"

samtools view \
    $home_dir/alignment/$annotate_file \
    -F 256 \
    -b \
    > $home_dir/alignment/$filtered_file


echo "## Getting properly paired reads"

proper_paired_file="${sample}_sorted_annotated_filtered_pp.bam"

samtools view \
    -bf 2x0 \
    $home_dir/alignment/$filtered_file \
    > $home_dir/alignment/$proper_paired_file

samtools index \
    $home_dir/alignment/$proper_paired_file

echo "## Running hmftools UMI deduplication"
#At the time of the creation of this pipeline the UMI deduplication tool
#from hmftools was called mark-dups. It has since been renamed redux
#which is available from:
#https://github.com/hartwigmedical/hmftools/releases/tag/redux-v2.0

hmftools_path="path/to/hmftools"
hmftools_unsorted_file="${sample}_hmftools_unsort.bam"

java \
    -jar $hmftools_path \
    -bam_file $home_dir/alignment/$proper_paired_file \
    -sample $sample \
    -ref_genome $reference_path/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -ref_genome_version V38 \
    -output_bam $home_dir/alignment/$hmftools_unsorted_file \
    -output_dir alignment/ \
    -write_stats \
    -umi_enabled \
    -umi_base_diff_stats 

echo "## Moving UMI stat files"

mv alignment/*tsv* umi_stats/

echo "## Sorting hmftools file"

hmftools_file="${sample}_hmftools.bam"

samtools sort \
    $home_dir/alignment/$hmftools_unsorted_file \
    > $home_dir/alignment/$hmftools_file

samtools index \
    $home_dir/alignment/$hmftools_file

echo "## Removing duplicates"

hmftools_rmdup_nosort_file="${sample}_rmdup.bam"

samtools view \
    $home_dir/alignment/$hmftools_file \
    -b \
    -F 1024 \
    > $home_dir/alignment/$hmftools_rmdup_nosort_file

samtools index \
    $home_dir/alignment/$hmftools_rmdup_nosort_file

echo "## Sorting BAM files according to mate"

hmftools_mate_file="${sample}_rmdup_mate.bam"

samtools sort \
    -n \
    $home_dir/alignment/$hmftools_rmdup_nosort_file \
    > $home_dir/alignment/$hmftools_mate_file

echo "## Creating BED file"

hmftools_bed_file="${sample}_rmdup.bed"

bedtools \
    bamtobed \
    -bedpe \
    -i $home_dir/alignment/$hmftools_mate_file \
    > $home_dir/bed/$hmftools_bed_file 

echo "## Creating 1bp BED files and calculating enrichment ratios"
#The R scripts are available at:
#https://github.com/CTrierMaansson/cfChIP_pembro/R/scripts/
R_script_paths="path/to/R_scripts"
#The regions defined as on_target and off_target regions of H3K4me3
#are available at:
#https://github.com/CTrierMaansson/cfChIP_pembro/reference/
target_path="path/to/reference/"

Rscript \
    --vanilla \
    $R_script_paths/Enrich_ratio_script.R \
    $home_dir/bed/$hmftools_bed_file \
    $target_path

echo "## Adding enrichment ratio to document"
hmftools_ratio_file="${sample}_rmdup_enrichment_ratio.txt"
hmftools_values="$(sed enrichment_ratios/$hmftools_ratio_file -n -e 2p | awk '{print $2}')"

if [[ ! -e $home_dir/enrichment_ratios/enrichment_ratios.txt ]]; then
    touch $home_dir/enrichment_ratios/enrichment_ratios.txt
    echo "ratios"  "sample" >> enrichment_ratios/enrichment_ratios.txt 
fi

echo $hmftools_values  $sample >> enrichment_ratios/enrichment_ratios.txt 

echo "## Creating BigWig files"

hmftools_bw_file="${sample}.bw"

hmftools_bigwig_file="${sample}.bigwig"

bamCoverage \
    --bam $home_dir/alignment/$hmftools_rmdup_nosort_file \
    -o $home_dir/bigwig/$hmftools_bw_file \
    --binSize 20 \
    --normalizeUsing CPM \
    --extendReads \
    --minMappingQuality 10
    
cp $home_dir/bigwig/$hmftools_bw_file $home_dir/bigwig/$hmftools_bigwig_file

echo "## Creating bedgraph files"

hmftools_bedgraph_file="${sample}.bedgraph"

bigWigToBedGraph \
    $home_dir/bigwig/$hmftools_bigwig_file \
    $home_dir/bedgraph/$hmftools_bedgraph_file
    
echo "## Getting fragment lengths"

Rscript 
    --vanilla 
    $R_script_paths/export_fragment_lengths.R \
    $home_dir/alignment/$hmftools_rmdup_nosort_file \
    $white_list

echo "# Preprocessing DONE"

