#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 32G
#SBATCH -c 1
#SBATCH --output=logs/preprocess%j.out
#SBATCH --time=36:00:00


sample=""
home_dir=""
fastq1=""
fastq2=""
reference_path=""
R_script_paths=""
genome_path=""
genome_prefix=""
fgbio_path=""
hmftools_path=""

usage() {
    echo "Usage: $0 Required arguments: -s SAMPLE_NAME -o OUTPUT_DIRECTORY -f FASTQ_R1_PATH -F FASTQ_R1_PATH -r REFERENCE_PATH -R R_PATH -g GENOME_PATH -p GENOME_PREFIX -b FGBIO_PATH -m HMFTOOLS_PATH"
    echo
    echo "Options:"
    echo "  -s SAMPLE   Sample name"
    echo "  -o DIRECTORY   Output directory"
    echo "  -f PATH   Path to uncompressed fastq R1"
    echo "  -F PATH   Path to uncompressed fastq R2"
    echo "  -r PATH   Path to reference data (https://github.com/CTrierMaansson/cfChIP_pembro/reference/)"
    echo "  -R PATH   Path to R scripts (https://github.com/CTrierMaansson/cfChIP_pembro/R/scripts/)"
    echo "  -g PATH   Path to bowtie indexed reference genome which also contain the genome FASTA and genome dictionary created with picard CreateSequenceDictionary"
    echo "  -p CHARACTER   Prefix of indexed genome, <prefix>.fa, from bowtie2-build"
    echo "  -b PATH   Path to fgbio .jar file (https://github.com/fulcrumgenomics/fgbio/releases)"
    echo "  -m PATH   Path to hmftools .jar file (https://github.com/hartwigmedical/hmftools/releases/tag/redux-v2.0)"
    echo "  -h        Show this help message"
}


while getopts ":s:o:f:F:r:R:g:p:b:m:h" opt; do
    case "$opt" in
        s)
            sample="$OPTARG"
            ;;
        o)
            home_dir="$OPTARG"
            ;;
        f)
            fastq1="$OPTARG"
            ;;
        F)
            fastq2="$OPTARG"
            ;;
        r)
            reference_path="$OPTARG"
            ;;
        R)
            R_script_paths="$OPTARG"
            ;;
        g)
            genome_path="$OPTARG"
            ;;
        p)
            genome_prefix="$OPTARG"
            ;;
        b)
            fgbio_path="$OPTARG"
            ;;
        m)
            hmftools_path="$OPTARG"
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
      -z "$fastq1" ||
      -z "$fastq2" ||
      -z "$reference_path" ||
      -z "$R_script_paths" ||
      -z "$genome_path" || 
      -z "$genome_prefix" ||
      -z "$fgbio_path" || 
      -z "$hmftools_path" ]]; then
    echo "Usage: $0 Required arguments: -s SAMPLE_NAME -o OUTPUT_DIRECTORY -f FASTQ_R1_PATH -F FASTQ_R1_PATH -r REFERENCE_PATH -R R_PATH -g GENOME_PATH -p GENOME_PREFIX -b FGBIO_PATH -m HMFTOOLS_PATH"
    exit 1
fi

# Reject compressed FASTQ files
if [[ "$fastq1" == *.gz ]]; then
    echo "Error: R1 input must not be gzipped: $fastq1" >&2
    exit 1
fi

if [[ "$fastq2" == *.gz ]]; then
    echo "Error: R2 input must not be gzipped: $fastq2" >&2
    exit 1
fi


require_file() {
    local path="$1"
    local option="$2"

    if [[ ! -f "$path" ]]; then
        echo "Error: $option must point to an existing regular file: $path" >&2
        exit 1
    fi
}

require_dir() {
    local path="$1"
    local option="$2"

    if [[ ! -d "$path" ]]; then
        echo "Error: $option must point to an existing directory: $path" >&2
        exit 1
    fi
}
fasta="${genome_prefix}.fa"

echo "Testing whether input files can be found"
echo "fastq_R1: ${fastq1}"
echo "fastq_R2: ${fastq2}"
echo "Reference genome FASTA: ${genome_path}/${fasta}"
echo "Reference whitelist: ${reference_path}/Hg38WhitelistV2.bed"
echo "Reference chromosome sizes: ${reference_path}/hg38_size.bed"
echo "H3K4me3 on-target file: ${reference_path}/H3K4me3_targets_reduced.bed"
echo "H3K4me3 off-target file: ${reference_path}/H3K4me3_offtargets_reduced.bed"
echo "Fragment length R script: ${R_script_paths}/export_fragment_lengths.R"
echo "Enrichment ratio R script: ${R_script_paths}/Enrich_ratio_script.R"
echo "TSS coverage R script: ${R_script_paths}/MANE_coverage.R"
echo "fgbio.jar file: ${fgbio_path}"
echo "hmftools.jar file: ${hmftools_path}"

require_file "$fastq1" "-f"
require_file "$fastq2" "-F"
require_file "$fgbio_path" "-b"
require_file "$hmftools_path" "-m"

require_dir "$home_dir" "-o"
require_dir "$reference_path" "-r"
require_dir "$R_script_paths" "-R"
require_dir "$genome_path" "-g"

if [[ ! -f "$genome_path/$fasta" ]]; then
        echo "Error: $fasta cannot be found in : $genome_path" >&2
        exit 1
fi

#The hg38 whitelist and other reference data 
#used for this manuscript is available at:
#https://github.com/CTrierMaansson/cfChIP_pembro/reference/

if [[ ! -f "$reference_path/Hg38WhitelistV2.bed" ]]; then
    echo "Error: Hg38WhitelistV2.bed cannot be found in: $reference_path" >&2
    exit 1
fi

if [[ ! -f "$reference_path/hg38_size.bed" ]]; then
    echo "Error: hg38_size.bed cannot be found in: $reference_path" >&2
    exit 1
fi

if [[ ! -f "$reference_path/H3K4me3_offtargets_reduced.bed" ]]; then
    echo "Error: H3K4me3_offtargets_reduced.bed cannot be found in: $reference_path" >&2
    exit 1
fi

if [[ ! -f "$reference_path/H3K4me3_targets_reduced.bed" ]]; then
    echo "Error: H3K4me3_targets_reduced.bed cannot be found in: $reference_path" >&2
    exit 1
fi

#The R scripts are available at:
#https://github.com/CTrierMaansson/cfChIP_pembro/R/scripts/

if [[ ! -f "$R_script_paths/export_fragment_lengths.R" ]]; then
    echo "Error: export_fragment_lengths.R cannot be found in: $R_script_paths" >&2
    exit 1
fi

if [[ ! -f "$R_script_paths/Enrich_ratio_script.R" ]]; then
    echo "Error: Enrich_ratio_script.R cannot be found in: $R_script_paths" >&2
    exit 1
fi

if [[ ! -f "$R_script_paths/MANE_coverage.R" ]]; then
    echo "Error: MANE_coverage.R cannot be found in: $R_script_paths" >&2
    exit 1
fi

echo "Test successful"

echo "# This is the preprocessing of the following sample: ${sample}"

echo "# Exporting results to ${home_dir}"

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
mkdir -p $home_dir/coverage
mkdir -p $home_dir/fragment_lengths


echo "## Running pre trimming fastQC"

fastqc \
    -o $home_dir/fastqc/ \
    -f fastq \
    $fastq1 \
    $fastq2

echo "## Running fastq trimming"

trim_galore \
    --fastqc \
    --fastqc_args "--outdir ${home_dir}/trimmed_fastqs/fastqc/" \
    --paired \
    -o $home_dir/trimmed_fastqs \
    $fastq1 \
    $fastq2

filename1="${fastq1##*/}"
sample1="${filename1%.fastq}"
filename2="${fastq2##*/}"
sample2="${filename2%.fastq}"

trimmed_R1_file="${sample1}_val_1.fq"
trimmed_R2_file="${sample2}_val_2.fq"
sam_file="${sample}.sam"

echo "## Aligning reads to ${genome_prefix}"

bowtie2 \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -x $genome_path/$genome_prefix \
    -p 4 \
    -1 $home_dir/trimmed_fastqs/$trimmed_R1_file \
    -2 $home_dir/trimmed_fastqs/$trimmed_R2_file \
    -S $home_dir/alignment/$sam_file 

bam_file="${sample}.bam"

echo "## Creating BAM file"

samtools view \
    -b \
    -L $reference_path/Hg38WhitelistV2.bed \
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

hmftools_unsorted_file="${sample}_hmftools_unsort.bam"

java \
    -jar $hmftools_path \
    -bam_file $home_dir/alignment/$proper_paired_file \
    -sample $sample \
    -ref_genome $genome_path/$fasta \
    -ref_genome_version V38 \
    -output_bam $home_dir/alignment/$hmftools_unsorted_file \
    -output_dir $home_dir/alignment/ \
    -write_stats \
    -umi_enabled \
    -umi_base_diff_stats 

echo "## Moving UMI stat files"

mv $home_dir/alignment/*tsv* $home_dir/umi_stats/

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


Rscript \
    --vanilla \
    $R_script_paths/Enrich_ratio_script.R \
    $home_dir/bed/$hmftools_bed_file \
    $reference_path/

echo "## Adding enrichment ratio to document"
hmftools_ratio_file="${sample}_rmdup_enrichment_ratio.txt"
hmftools_values="$(sed $home_dir/enrichment_ratios/$hmftools_ratio_file -n -e 2p | awk '{print $2}')"

if [[ ! -e $home_dir/enrichment_ratios/enrichment_ratios.txt ]]; then
    touch $home_dir/enrichment_ratios/enrichment_ratios.txt
    echo "ratios"  "sample" >> $home_dir/enrichment_ratios/enrichment_ratios.txt 
fi

echo $hmftools_values  $sample >> $home_dir/enrichment_ratios/enrichment_ratios.txt 

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
    
echo "## Getting TSS coverage across genes"


Rscript \
    --vanilla 
    $R_script_paths/MANE_coverage.R \
    $home_dir/bedgraph/$hmftools_bedgraph_file 
    
echo "## Getting fragment lengths"


Rscript \
    --vanilla 
    $R_script_paths/export_fragment_lengths.R \
    $home_dir/alignment/$hmftools_rmdup_nosort_file \
    $reference_path/hg38_size.bed

echo "# Preprocessing DONE"

