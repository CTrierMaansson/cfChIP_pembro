args = commandArgs(trailingOnly=TRUE)

bedgraph_sample <- args[1]

bname <- basename(bedgraph_sample)
patient_short <- gsub(".bedgraph","",bname)


library(dplyr)
library(GenomicFeatures)
library(GenomicRanges)
print("Defining MANE TSS regions")
MANE_genes <- GenomicFeatures::genes(TxDb.Hsapiens.NCBI.hg38.MANE::TxDb.Hsapiens.NCBI.hg38.MANE)
MANE_genes_df <- as.data.frame(MANE_genes) %>% 
  mutate(TSS = ifelse(strand == "-", end, start)) %>% 
  dplyr::select(-c(start, end, width)) %>% 
  mutate(start = TSS-3000) %>% 
  mutate(end = TSS+3000)  %>% 
  mutate(GENEID = unlist(lapply(stringr::str_split(gene_id, "[.]"), "[[", 1)))

print("Adding SYMBOLs")
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

MANE_genes_df <- MANE_genes_df %>% 
  left_join(gene_id_df, by = "GENEID") %>% 
  dplyr::filter(!is.na(SYMBOL))


MANE_TSS_regions <- with(MANE_genes_df,
                         GRanges(seqnames,
                                 IRanges(start,end),
                                 strand,
                                 TSS = TSS,
                                 geneId = gene_id,
                                 geneId_short = GENEID,
                                 SYMBOL = SYMBOL))


track_df <- function(x,y,sample){
  library(plyr)
  library(dplyr)
  library(plyranges)
  print("Reading bedGraph")
  patient_gr <- read_bed_graph(x)
  print("Subsetting bedGraph")
  patient_gr_subset <- subsetByOverlaps(patient_gr,y)
  print("Analysing windows")
  df <- data.frame(pos = seq(-3000,3000,20))
  for(i in 1:length(y)){
    if(i%%100 == 0){
      print(paste("Progress is",round(i/length(y)*100,1),"%"))
    }
    gr_sub <- y[i]
    bg_sub <- suppressWarnings(subsetByOverlaps(patient_gr_subset,gr_sub))
    TSS <- gr_sub$TSS
    relative_pos <- round_any(((end(bg_sub)-start(bg_sub))/2)+start(bg_sub)-TSS,20)
    st <- as.character(strand(gr_sub))
    if(st == "-"){
      relative_pos <- rev(relative_pos)
    }
    ddf <- data.frame(pos = relative_pos,
                      score = elementMetadata(bg_sub)$score) %>% 
      dplyr::filter(pos <= 3000 & pos >= -3000)
    colnames(ddf) <- c("pos", elementMetadata(gr_sub)$SYMBOL)
    df <- full_join(df,ddf, by = "pos")
    df[,i+1][is.na(df[,i+1])] <- 0
  } 
  print(paste("Progress is",round(i/length(y)*100,1),"%"))
  df$sample <- paste(sample)
  return(df)
}
print("Extracting TSS coverage")
patient_df <- track_df(bedgraph_sample,
                       MANE_TSS_regions,
                       patient_short)

root_path <- gsub("bedgraph",
                  "",
                  dirname(bedgraph_sample))

print("Writting TSS profile file")
write.table(patient_df, 
            file = paste0(root_path,"coverage/",patient_short,"_H3K4me3_MANE_TSS_coverage.txt"), 
            sep = "\t",
            col.names = T,
            row.names = F)
