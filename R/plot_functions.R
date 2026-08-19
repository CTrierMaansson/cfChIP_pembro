# Functions used for plot generation

#### Loading libraries####
library(dplyr)
library(ggpubr)
library(ggplot2)
library(viridisLite)


####Defining plot theme####
th <- theme(
    legend.background = element_rect(),
    plot.title = element_text(angle = 0,
                              size = 12,
                              face = "bold",
                              vjust = 1,
                              hjust = 0.5),
    plot.caption = element_text(angle = 0,
                                size = 10,
                                vjust = 1,
                                hjust = 0.37),
    axis.text.x = element_text(angle = 0, 
                               size = 10),
    axis.text.y = element_text(angle = 0,
                               size = 10),
    axis.title = element_text(size = 10,
                              face = "bold"),
    axis.title.x = element_text(size = 10,
                                face = "bold"),
    axis.title.y = element_text(size = 10,
                                face = "bold"),
    axis.line = element_line(color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10,
                                face = "bold"))

####Defining functions####

# H3K4me3 enrichment vs. TF
CRE_ctDNA_single <- function(counts, #Gene level cfChIP enrichment returned by <functionname>
                             meta, #cfChIP metadata for each sample
                             SYMBOL, #Gene SYMBOL to plot enrichment for
                             er = NULL, #Enrichment ratio cutoff
                             lib_size = NULL #Number of cfChIP fragments cutoff
                             ) {
    #Filtering samples
    meta <- meta %>% 
        filter(study_name == "Lung cancer")
    meta <- filter_samples(meta = meta,
                           er = er,
                           lib_size = lib_size)
    meta_sele <- meta %>% 
        dplyr::select(sample, ctDNA,sample_type) %>% 
        mutate(sample_type = ifelse(sample_type == "Baseline","BL","Tx"))
    #Selecting relevant samples
    counts <- counts[SYMBOL,meta$sample]
    df <- as.data.frame(counts)
    df$sample <- rownames(df)
    rownames(df) <- NULL
    df <- df %>% 
        inner_join(meta_sele, 
                   by ="sample")
    #Performing correlation analysis
    res <- cor.test(df$counts, 
                    df$ctDNA,
                    method = "spearman")
    if(res$estimate < 0){
        x_r <- 10
    }else{
        x_r <- 1
    }
    #Formatting Spearman's r and p-value
    spear <- sprintf("%.2f",res$estimate)
    ps <- format_p(res$p.value)
    #Creating plot
    gg <- ggplot(df, 
                 aes(x = ctDNA,
                     y = counts))+
        geom_point(aes(shape = sample_type),
                   size = 3, color = "#751F58")+
        scale_y_continuous(trans = "log10",
                           expand = c(0,0.1))+
        scale_x_continuous(trans = "log10",
                           expand = c(0,0.1))+
        labs(x = "Tumor fraction [%]",
             y = "H3K4me3 enrichment",
             title = SYMBOL)+
        spear_p(spear = spear,
                p = ps,
                y = Inf,
                x = x_r,
                h = 0,
                v = 1)+
        geom_smooth(method = "lm", se = T,
                    color = "black")+
        scale_shape_discrete("Sample type   ")+
        theme_classic()+
        th+
        theme(legend.position = "top")
    return(gg)
}

# GSEA results as bar graph
GSEA_barplot <- function(gsea_res, #GSEA results returned by <functionname>
                         n_terms, #Number of terms to display in plot
                         title = NULL, #Add title to plot
                         color = "#751F58" #Color of bars
                         ){
    df <- gsea_res
    #Selecting relevant gene sets
    df <- df %>% 
        arrange(desc(abs(NES))) %>% 
        slice_head(n = n_terms)
    #Simplifying gene set names
    df$pathway <- gsub("_"," ",df$pathway)
    df$pathway <- gsub("HALLMARK ","",df$pathway)
    df <- df %>% 
        mutate(pathway = factor(pathway, levels = rev(df$pathway)))
    gg <- ggplot(df, 
                 aes(x = pathway, 
                     y = NES))+
        geom_point(size = 4, 
                   color = color)+
        geom_col(stat = "identity",
                 width = 0.1,
                 color = color,
                 fill = color)+
        theme_classic()+
        geom_hline(yintercept = 0,
                   colour = "black",
                   linewidth = 0.8,
                   linetype = "solid")+
        labs(y = "nES",
             x = "")+
        coord_flip()+
        th
    #Adding title
    if(!is.null(title)){
        gg <- gg+
            labs(title = title)
    }
    return(gg)
}

#PCA plot
basic_PCA <- function(counts, #Gene level cfChIP enrichment returned by <functionname>
                      meta, #cfChIP metadata for each sample
                      genes, #Gene SYMBOLS which the PCA is based upon
                      samples, #List of samples/patients, where the names of the 
                               #list entries refer to the sample groups 
                               #determining the color in the plot
                      er = NULL, #Enrichment ratio cutoff
                      lib_size = NULL, #Number of cfChIP fragments cutoff
                      color = NULL, #Definition of colors to use, otherwise viridis(begin = 0.3, end = 0.9) is used
                      title = NULL #Add title to plot
                      ){
    # Selecting samples to plot
    all_samples <- unlist(samples)
    meta <- meta %>% 
        filter(study_name != "Replicate") %>% 
        filter(sample %in% all_samples | patient %in% all_samples)
    #Filtering meta information
    meta <- filter_samples(meta = meta,
                           er = er,
                           lib_size = lib_size)
    if(any(grepl("_",all_samples))){
        counts_sele <- counts[genes,meta$sample]
    }else{
        counts_sele <- counts[genes,meta$patient]
    }
    #Performing PCA
    counts_sele <- counts_sele %>% t() %>% unique()
    counts_sele <- counts_sele[,colSums(counts_sele) != 0]
    pca_res <- prcomp(x = counts_sele, 
                      scale = TRUE,
                      center = TRUE)
    sum_res <- summary(pca_res)
    #Formatting results for the plot
    PC1 <- paste0("PC-1 (",
                  round(sum_res$importance[2]*100,
                        digits = 2),
                  "%)")
    PC2 <- paste0("PC-2 (",
                  round(sum_res$importance[5]*100,
                        digits = 2),
                  "%)")
    pca_res_df <- as.data.frame(pca_res$x) %>% 
        dplyr::select(PC1, PC2)
    pca_res_df$sample <- rownames(pca_res_df)
    #Merging PCA results and samples
    sample_df <- stack(samples)
    colnames(sample_df) <- c("sample","group")
    df <- pca_res_df %>% 
        left_join(sample_df, by = "sample")
    #Defining colors
    if(is.null(color)){
        vals <- viridis(begin = 0.3, end = 0.9,n = length(samples))
        names(vals) <- names(samples)
    }else{
        vals <- color
        names(vals) <- names(samples)
    }
    #Creating plot
    gg <- ggplot(data = df,
                 aes(x = PC1, 
                     y = PC2, 
                     color = group))+
        geom_point(size = 3.5)+
        labs(x = PC1, 
             y = PC2)+
        scale_color_manual(name = "Sample   ",
                           values = vals)+
        theme_classic()+
        th+
        theme(legend.position = "top")
    #Adding title
    if(!is.null(title)){
        gg <- gg+
            labs(title = title)
    }
    return(gg)
}

#Plotting TSS enrichment boxplots for selected genes
TSS_counts_selected_genes <- function(counts, #Gene level cfChIP enrichment returned by <functionname>
                                      meta, #cfChIP metadata for each sample
                                      genes, #Gene SYMBOLS which should be plotted.
                                      samples, #List of samples, where the names of the 
                                      #list entries refer to the sample groups 
                                      #determining the x-value in the boxplot
                                      er = NULL, #Enrichment ratio cutoff
                                      lib_size = NULL, #Number of cfChIP fragments cutoff
                                      color = NULL, #Definition of colors to use, otherwise viridis(begin = 0.3, end = 0.9) is used
                                      concat = F, #If multiple genes are listed,
                                                  #concat = T will plot the average enrichment,
                                                  #if concat = F, indivual plots will be created for each gene 
                                      title = NULL #Add plot title
                                      ){
    #Filtering genes and samples
    g_df <- data.frame(sig = genes) %>% 
        rename_AC(char = "sig")
    all_samples <- unlist(samples)
    meta <- meta %>% 
        filter(study_name != "Replicate") %>% 
        filter(sample %in% all_samples)
    meta <- filter_samples(meta = meta,
                           er = er,
                           lib_size = lib_size)
    counts_sele <- counts[genes,meta$sample]
    df <- as.data.frame(counts_sele)
    #Formatting enrichment and gene information 
    if(length(genes) == 1){
        df$SYMBOL <- genes
        df$sample <- rownames(df)
        rownames(df) <- NULL
        colnames(df)[1] <- "signal"
        df <- df %>% 
            dplyr::select(sample,signal,SYMBOL) %>% 
            rename_AC(char = "SYMBOL")
    }else{
        df$SYMBOL <- rownames(df)
        rownames(df) <- NULL
        df <- df %>% 
            pivot_longer(-SYMBOL) %>% 
            dplyr::rename("sample" = "name") %>% 
            dplyr::rename("signal" = "value") %>% 
            dplyr::select(sample,signal,SYMBOL) %>% 
            rename_AC(char = "SYMBOL") %>% 
            mutate(SYMBOL = factor(SYMBOL, levels = g_df$sig))
    }
    #Calculating average enrichment if concat = TRUE
    if(concat){
        mean_vals <- counts_sele %>% 
            t() %>% 
            as.data.frame() %>% 
            rowMeans()
        df <- data.frame(sample = names(mean_vals),
                         signal = mean_vals)
    }
    #Comparing enrichment between groups, defined by "samples"
    sample_df <- stack(samples)
    colnames(sample_df) <- c("sample","group")
    df <- df %>% 
        left_join(sample_df, by = "sample")
    compares <- combn(names(samples),m = 2,
                      simplify = F)
    if(concat){
        stat_test <- df %>% 
            t_test(signal ~ group,
                   comparisons = compares,
                   paired = FALSE) %>% 
            add_xy_position(x = "group",
                            data = df,
                            formula = signal ~ group,
                            scales = "free_y",
                            step.increase = 0.12) 
    }else{
        stat_test <- df %>% 
            group_by(SYMBOL) %>% 
            t_test(signal ~ group,
                   comparisons = compares,
                   paired = FALSE) %>% 
            add_xy_position(x = "group",
                            data = df,
                            formula = signal ~ group,
                            scales = "free_y",
                            step.increase = 0.12) 
    }
    #Defining colors
    if(is.null(color)){
        vals <- viridis(begin = 0.3, end = 0.9,n = length(samples))
        names(vals) <- names(samples)
    }else{
        vals <- color
        names(vals) <- names(samples)
    }
    #Creating plot for gene enrichment average
    if(concat){
        gg <- ggplot(df,
                     aes(x = group,
                         y = signal,
                         fill = group))+
            geom_boxplot(outlier.alpha = 0)+
            scale_fill_manual(values = vals)+
            theme_classic()+
            labs(y = "Mean H3K4me3 enrichment",
                 x = "")+
            scale_y_continuous(expand = c(0.15,0))+
            g_brac(stat_test)+
            th+
            theme(legend.position = "none")
    }
    #Creating plot for individual genes
    else{
        gg <- ggplot(df,
                     aes(x = group,
                         y = signal,
                         fill = group))+
            geom_boxplot(outlier.alpha = 0)+
            scale_fill_manual(values = vals)+
            theme_classic()+
            labs(y = "H3K4me3 enrichment",
                 x = "")+
            facet_wrap(~SYMBOL, scales = "free",nrow = 1)+
            scale_y_continuous(expand = c(0.15,0))+
            g_brac(stat_test)+
            th+
            theme(legend.position = "none",
                  strip.background = element_blank(),
                  strip.text = element_text(face = "bold",
                                            size = 12),
                  panel.background = element_rect(fill = "white",
                                                  color = "white"),
                  plot.background = element_rect(fill = "white",
                                                 color = "white"))
    }
    #Adding title
    if(!is.null(title)){
        gg <- gg +
            labs(title = title)
    }
    return(gg)
}