args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(GenomicRanges)

bam_file <- args[1]
size_file <- args[2]
bname <- basename(bam_file)
sample <- gsub("_rmdup.bam","",bname) #Defining sample name

message(paste0("Reading peak file for sample: ",sample))

root_path <- dirname(dirname(bam_file))

sample_repo <- paste0(root_path,"/macs/",sample,"/")
peak_file <- paste0(sample_repo,sample,"_peaks.broadPeak")
peaks <- fread(peak_file)
message("Creating GRanges")
peaks_gr <- with(peaks, 
                 GRanges(V1,
                         IRanges(start = V2, 
                                 end = V3)))
extended_peaks <- peaks_gr+1000 #Adding 1kb to each side of the peaks to ensure peak regions are not classified as background
reduced_peaks <- reduce(extended_peaks)

message("Reading hg38 chromosome sizes")
hg38 <- fread(size_file)
hg38_gr <- with(hg38,
                GRanges(V1, 
                        IRanges(start = 1, 
                                end = V2)))

message("Finding regions not covered by peaks")
gaps_gr <- GenomicRanges::setdiff(hg38_gr,reduced_peaks)

gr_to_bed <- function(x,y, z = NULL){
    df <- data.frame(chr = seqnames(x),
                     start = start(x),
                     end = end(x),
                     strand = strand(x))
    if(!is.null(z)){
        for (i in 1:length(z)){
            df[,(4+i)] <- unlist(elementMetadata(x)[,z[i]])
            colnames(df)[(4+i)] <- z[i]
        }
    }
    write.table(df, 
                y, 
                sep="\t", 
                col.names=F,
                row.names = FALSE, 
                append = F, 
                quote = FALSE)
}

region_bed_path <- paste0(root_path,"/IchorCNA/background_regions/")

message("Writing BED files")
gr_to_bed(x = gaps_gr ,
          y = paste0(region_bed_path,sample,"_background_regions.bed"))
message("DONE")