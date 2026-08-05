args = commandArgs(trailingOnly=TRUE)

bam_sample <- args[1]

white_list_path <- args[2]

library(Rsamtools)
library(GenomicAlignments)
library(dplyr)
library(data.table)

#Defining the sizes of hg38 as well as relevant regions
hg38_gr <- with(fread(white_list_path),
                GRanges(V1,
                        IRanges(V2,V3)))
param <- ScanBamParam(which = hg38_gr,
                      what = c("isize","pos"))

print("Reading BAM file")
print(paste(bam_sample))
bam <- readGAlignments(file = bam_sample,
                       param = param)
bam <- bam[!is.na(mcols(bam)$isize)] #Filtering NA fragment lengths

print("Creating fragment length table")
bam_single <- bam[mcols(bam)$isize > 0] #Ensures we only count fragments once
df_length <- as.data.frame(table(mcols(bam_single)$isize)) %>% 
  mutate(prop = Freq/sum(Freq))

print("Creating file name")
bname <- basename(bam_sample)
patient_sample <- gsub("_rmdup.bam","",bname)
print("Sample name:")
print(paste(patient_sample))

print("Crating data.frame of fragment lengths")
colnames(df_length) <- c("fragment_length","count","fraction")
df_length <- df_length %>% 
  mutate(sample = patient_sample) %>% 
    mutate(fragment_length = as.numeric(fragment_length)) %>% 
    filter(fragment_length < 501)

root_path <- gsub("alignment/",
                  "",
                  gsub(bname,
                       "",
                       bam_sample))

print("Writting fragment length file")
write.table(df_length, 
            file = paste0(root_path,"fragment_lengths/",patient_sample,".txt"), 
            sep = "\t",
            col.names = T,
            row.names = F)