args = commandArgs(trailingOnly=TRUE)
library(dplyr)

peak_file <- args[1]

load_peak <- function(x){
    library(data.table)
    library(GenomicAlignments)
    peak.file <- fread(x)
    peak.gr <- with(peak.file, 
                    GRanges(V1, 
                            IRanges(start = V2, 
                                    end = V3)))
    return(peak.gr)
    
}

ID_SYMBOL <- function(){
    mane_TxDb <- TxDb.Hsapiens.NCBI.hg38.MANE::TxDb.Hsapiens.NCBI.hg38.MANE
    mane_genes <- GenomicFeatures::genes(mane_TxDb)
    gene_ids <- unlist(lapply(strsplit(mane_genes$gene_id,".",fixed = TRUE),"[[",1))
    gene_id_df <- ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                    keys = gene_ids,
                                    keytype = "GENEID",
                                    columns = c("SYMBOL"))
    return(gene_id_df)
}
gene_id_df <- ID_SYMBOL()

bname <- basename(peak_file)
sample <- gsub("_peaks.broadPeak","",bname) #Defining sample name

#Defining output path
macs_path <- gsub(sample,
                  "",
                  dirname(peak_file))
output_path <- paste0(macs_path,"annotated/")

annotated_file <- paste0(output_path, sample, "_MANE_peaks.txt") #Defining annotated file name
stats_file <- paste0(output_path, sample, "_MANE_stats.txt") #Defining peak stats file name

print("Loading peaks")
peaks <- load_peak(peak_file) #Loading the peaks
print("Annotating peaks")
mane_TxDb <- TxDb.Hsapiens.NCBI.hg38.MANE::TxDb.Hsapiens.NCBI.hg38.MANE
annotated_peaks <- ChIPseeker::annotatePeak(peaks, #Annotating the peaks
                                            TxDb = mane_TxDb,
                                            level = "transcript",
                                            annoDb = "org.Hs.eg.db")


print("Calculating peak stats")
stats_df <- ChIPseeker::getAnnoStat(annotated_peaks) #Getting the annotation stats
stats_df <- apply(stats_df,2,as.character)

print("Adding SYMBOL information")
peak_df <- as.data.frame(annotated_peaks)
gene_ids <- unlist(lapply(strsplit(peak_df$geneId,".",fixed = TRUE),"[[",1))
peak_df$GENEID <- gene_ids

peak_df <- peak_df %>% 
    dplyr::select(-c(geneId,SYMBOL)) %>% 
    left_join(gene_id_df, by ="GENEID") %>% 
    filter(!is.na(SYMBOL))

print("Writting peak annotation stats")
write.table(
    stats_df,
    stats_file,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE)
print("Writting annotated peaks")
peak_df <- apply(peak_df,2,as.character)
write.table(
    peak_df,
    annotated_file,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE)