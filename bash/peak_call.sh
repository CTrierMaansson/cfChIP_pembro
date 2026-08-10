#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH -c 1
#SBATCH --output=logs/peak_call%j.out
#SBATCH --time=08:00:00


sample=""
home_dir=""
R_script_paths=""

usage() {
    echo "Usage: $0 Required arguments: -s SAMPLE_NAME -o OUTPUT_DIRECTORY -R R_PATH"
    echo
    echo "Options:"
    echo "  -s SAMPLE   Sample name matching a BAM file alignment/<SAMPLE>_rmdup.bam"
    echo "  -o DIRECTORY   Output directory"
    echo "  -R PATH   Path to R scripts (https://github.com/CTrierMaansson/cfChIP_pembro/R/scripts/)"
}

while getopts ":s:o:R:h" opt; do
    case "$opt" in
        s)
            sample="$OPTARG"
            ;;
        o)
            home_dir="$OPTARG"
            ;;
        R)
            R_script_paths="$OPTARG"
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
      -z "$R_script_paths" ]]; then
    echo "Usage: $0 Required arguments: -s SAMPLE_NAME -o OUTPUT_DIRECTORY -R R_PATH"
    exit 1
fi


require_dir() {
    local path="$1"
    local option="$2"

    if [[ ! -d "$path" ]]; then
        echo "Error: $option must point to an existing directory: $path" >&2
        exit 1
    fi
}
cfChIP_file="${sample}_rmdup.bam"

echo "Testing whether input files can be found"
echo "cfChIP BAM file: $home_dir/alignment/$cfChIP_file"
echo "Annotate peak R script: ${R_script_paths}/annotate_peaks.R"

require_dir "$home_dir" "-o"
require_dir "$R_script_paths" "-R"

if [[ ! -f "$R_script_paths/annotate_peaks.R" ]]; then
    echo "Error: annotate_peaks.R cannot be found in: $R_script_paths" >&2
    exit 1
fi

if [[ ! -f "$home_dir/alignment/$cfChIP_file" ]]; then
    echo "Error: cannot find $cfChIP_file in $home_dir/alignment/" >&2
    exit 1
fi

echo "Test successful"

echo "# This is the peak calling and annotation of the following sample: ${sample}"

eval "$(conda shell.bash hook)"
conda activate cfChIP_pembro

mkdir -p $home_dir/macs
mkdir -p $home_dir/macs/$sample
mkdir -p $home_dir/macs/annotated

echo "# Calling peaks"

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

echo "# Annotating peaks"

#The R scripts are available at:
#https://github.com/CTrierMaansson/cfChIP_pembro/R/scripts/

Rscript \
    --vanilla \
    $R_script_paths/annotate_peaks.R \
    $home_dir/macs/$sample/$peak_file
