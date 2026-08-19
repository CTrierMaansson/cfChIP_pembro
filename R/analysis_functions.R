# Functions used for data analysis


#### Loading libraries####
library(dplyr)

####Defining functions####

#Group patients based on cfChIP enrichment
dichotomize_patients_enrichment <- function(counts, #Gene level cfChIP enrichment returned by <functionname>
                                            meta, #cfChIP metadata for each sample
                                            SYMBOL, #Gene SYMBOL where the enrichment will separate patients
                                            samples, #samples which should be used for the analysis
                                            er = NULL, #Enrichment ratio cutoff
                                            lib_size = NULL #Number of cfChIP fragments cutoff
                                            ){
    #Filter samples
    meta <- meta %>% 
        dplyr::filter(study_name == "Lung cancer") %>% 
        dplyr::filter(sample %in% samples)
    meta <- filter_samples(meta = meta,
                           er = er,
                           lib_size = lib_size)
    #Select genes and samples
    counts_sele <- counts[SYMBOL,meta$sample]
    df <- as.data.frame(counts_sele)
    df$sample <- rownames(df)
    rownames(df) <- NULL
    colnames(df)[1] <- "signal"
    meta_sele <- meta %>% 
        dplyr::select(sample, patient)
    df <- df %>% 
        left_join(meta_sele, by = "sample")
    #Define median and group samples
    med_signal <- median(df$signal)
    df <- df %>% 
        mutate(divider = ifelse(signal >= med_signal,1,0))
    return(df)
}

