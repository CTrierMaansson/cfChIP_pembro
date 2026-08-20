# Various functions used for data formatting or variable definitions

#### Loading libraries####
library(org.Hs.eg.db)
library(TxDb.Hsapiens.NCBI.hg38.MANE)
library(GenomicFeatures)
library(GenomicAlignments)
library(ensembldb)
library(dplyr)
library(ggplot2)
library(ggpubr)

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

#Formatting p-value to use in plots
format_p <- function(x,
                     cutoff=0.001,
                     digits=2){
    x <- signif(x,digits)
    vapply(x,function(z){
        if(is.na(z)){
            NA_character_
        } else if(z<= 2.2e-16){
            
            "< 2.2%*% 10^-16"
            
        }else if(abs(z) < cutoff){
            s <- formatC(z,
                         format = "e",
                         digits = digits-1,
                         drop0trailing = TRUE)
            parts <- strsplit(s,"e",fixed = TRUE)[[1]]
            coefficient <- parts[1]
            exponent  <- as.integer(parts[2])
            paste0(coefficient," %*% 10^",exponent)
        }else{
            formatC(z,
                    format = "fg",
                    digits = digits,
                    drop0trailing = TRUE)
        }
    },
    character(1))
}

#Adding p-value to plots comparing two groups
g_brac <- function(stat_test){
    stat_test <- stat_test %>% 
        mutate(p_label = ifelse(p<2.2e-16,
                                "p < 2.2%*% 10^-16",
                                paste0("p == ",format_p(p,
                                                        cutoff = 0.001,
                                                        digits = 2))))
    geom_bracket(data = stat_test,
                 aes(xmin = group1,
                     xmax = group2,
                     y.position = y.position,
                     label = p_label),
                 type = "expression",
                 tip.length = 0.02,
                 bracket.nudge.y = 0.1,
                 inherit.aes = FALSE)
}

#Formatting p-values and Spearman's r in plots
spear_p <- function(spear,p,y,x,h,v,nudge_x = 0){
    if(p == "< 2.2%*% 10^-16"){
        p_lab <- "p < 2.2%*% 10^-16"
    }else{
        p_lab = paste0("p == ",
                       p)
    }
    list(annotate("text",
                  y = y,
                  x = x,
                  label = paste0("r == '",
                                 spear,"'"),
                  parse = TRUE,
                  size = 4,
                  hjust = h,
                  vjust = v),
         annotate("text",
                  y = y,
                  x = x,
                  label = p_lab,
                  parse = TRUE,
                  size = 4,
                  hjust = h,
                  vjust = v+1.2))
}
#Defining peak gene annotations
new_annotation <- function(x){
    res <- unlist(lapply(x, FUN = function(annotation){
        if(grepl("Exon", annotation)){
            if(grepl("exon 1", annotation)){
                return("1st Exon")
            }
            else{
                return("Other Exon")
            }
        }
        if(grepl("Intron",annotation)){
            if(grepl("intron 1",annotation)){
                return("1st Intron")
            }
            else{
                return("Other Intron")
            }
        }
        return(annotation)
    }))
    return(res)
}