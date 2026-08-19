# Various functions used for data formatting or variable definitions

#### Loading libraries####
library(org.Hs.eg.db)
library(TxDb.Hsapiens.NCBI.hg38.MANE)
library(GenomicFeatures)
library(GenomicAlignments)
library(ensembldb)
library(dplyr)

####Defining functions####

# Collecting MANE transcript and Gene SYMBOL info
ID_SYMBOL <- function(){
    #Get MANE transcripts (https://github.com/CTrierMaansson/TxDb.Hsapiens.NCBI.hg38.MANE)
    mane_TxDb <- TxDb.Hsapiens.NCBI.hg38.MANE::TxDb.Hsapiens.NCBI.hg38.MANE
    mane_genes <- GenomicFeatures::genes(mane_TxDb)
    gene_ids <- unlist(lapply(strsplit(mane_genes$gene_id,".",fixed = TRUE),"[[",1))
    #Merging SYMBOL and GENEID
    gene_id_df <- ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                    keys = gene_ids,
                                    keytype = "GENEID",
                                    columns = c("SYMBOL"))
    return(gene_id_df)
}
gene_id_df <- ID_SYMBOL()

# Filtering samples based on meta information
filter_samples <- function(meta,
                           er = NULL,
                           lib_size = NULL){
    if(!is.null(er)){
        meta <- meta %>% 
            dplyr::filter(enrichment > er)
    }
    if(!is.null(lib_size)){
        meta <- meta %>% 
            dplyr::filter(fragments > lib_size)
    }
    return(meta)
}

#Renaming AC013461.1 to MAP3K20
#Because MAP3K20 is defined as AC013461.1 in ensembldb we used this function
# to manually rename the gene in plots
rename_AC <- function(df,char){
    idx <- which(colnames(df) == char)
    colnames(df)[idx] <- "sig"
    df <- df %>% 
        mutate(sig = ifelse(sig == "AC013461.1","MAP3K20",sig))
    colnames(df)[idx] <- char
    return(df)
}
