args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
    stop(
        "Usage: Rscript Enrich_ratio_script.R <bed_file> <reference_path>",
        call. = FALSE
    )
}

`%ni%` <- Negate(`%in%`)

bed_sample <- args[1]

target_path <- args[2]

bname <- basename(bed_sample)

patient_sample <- unlist(lapply(stringr::str_split(bname, "[.]"), "[[", 1))
file_type <- unlist(lapply(stringr::str_split(bname, "[.]"), "[[", 2))

if (file_type %ni% c("bed")) {
    stop("File should be .bed", call.=FALSE)
}

enrichment_ratio_file <- paste0(patient_sample, "_enrichment_ratio.txt")
single_bp_file <- paste0(patient_sample, "_1bp.bed")

root_path <- dirname(dirname(bed_sample))

#x name of .bed file of off-target or on-target bins returned by off_target_gr() or reduce_target_table()
import_target <- function(x){
    library(dplyr)
    library(GenomicRanges)
    df <- as.data.frame(
        read.table(x,
                   header = FALSE,
                   sep="\t",
                   stringsAsFactors=FALSE,
                   quote="")
    )
    colnames(df) <- c("chr","start","end","range")
    gr <- GRanges(seqnames = df$chr, 
                  ranges = IRanges(start = df$start,
                                   end = df$end),
                  strand = "*")
    return(gr)
}

regular_chr <- function(){
    ends <- c(1:22,"X","Y")
    chr_ends <- paste0("chr",ends)
    return(chr_ends)
}

#x name of BED file representing reads from a cfChIP file
read_cfChIP_bed <- function(x){
    library(GenomicRanges)
    library(dplyr)
    chrs <- regular_chr()
    bed <- as.data.frame(
        read.table(x,
                   header = FALSE,
                   sep="\t",
                   stringsAsFactors=FALSE,
                   quote="")) %>% 
        filter(V1 == V4) %>% #Ensuring the read pairs are aligned to the same chromosome
        filter(V1 %in% chrs) #Only using reads aligned to the regular chromosomes
    gr <- GRanges(seqnames = bed$V1,
                  ranges = IRanges(start = bed$V2, #Start of cfChIP fragment
                                   end = bed$V6), #End of cfChIP fragment
                  strand = "*")
    return(gr)
}

#x bedfile Granges returned by read_cfChIP_bed()
short_bed <- function(x){
    library(GenomicRanges)
    starts <- start(x)
    ends <- end(x)
    widths <- width(x)
    poss <- starts+round(widths/2) #Middle of fragment
    gr <- GRanges(seqnames = seqnames(x),
                  ranges = IRanges(start = poss,
                                   end = poss),
                  strand = "*")
    return(gr)
}

#x GRanges object of 1bp reduced reads from a H3K4me3 cfChIP file returned by short_bed()
#y GRanges object of on-target bins returned by import_target()
#z GRanges object of off-target bins returned by import_target()
on_off_ratio <- function(x,y,z){
    library(IRanges)
    target_overlaps <- sum(countOverlaps(y,x))
    offtarget_overlaps <- sum(countOverlaps(z,x))
    norm_target_overlaps <- target_overlaps/(sum(width(y)))
    norm_offtarget_overlaps <- offtarget_overlaps/(sum(width(z)))
    ratio <- norm_target_overlaps/norm_offtarget_overlaps
    df <- data.frame(normalized_ratio = ratio,
                     ratio = target_overlaps/offtarget_overlaps)
    return(df)
}

print("Importing on-target .bed file")

target_gr <- import_target(paste0(target_path,
                                  "H3K4me3_targets_reduced.bed"))

print("Importing off-target .bed file")

offtarget_gr <- import_target(paste0(target_path,
                                     "H3K4me3_offtargets_reduced.bed"))

print("Importing cfChIP .bed file")

sample_bed <- read_cfChIP_bed(bed_sample)

print("Creating 1bp cfChIP GRanges object")

sample_short_bed <- short_bed(sample_bed)

single_bp <- data.frame(chr = seqnames(sample_short_bed),
                        start = start(sample_short_bed),
                        end = end(sample_short_bed))

print("Writing cfChIP 1bp .bed file")

write.table(single_bp, 
            file = paste0(root_path,
                          "/reduced_fragments/",
                          single_bp_file), 
            sep = "\t",
            col.names = T)


res <- on_off_ratio(sample_short_bed, target_gr,offtarget_gr)

print("Writing enrichment ratio file")

write.table(res, 
            file = paste0(root_path,
                          "/enrichment_ratios/",
                          enrichment_ratio_file), 
            sep = "\t",
            col.names = T)

print("Enrichment ratio analysis done")