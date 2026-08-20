# Functions used for plot generation

#### Loading libraries####
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(viridisLite)
library(surivival)
library(survminer)
library(reshape2)
library(rstatix)


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
CRE_ctDNA_single <- function(counts, #Gene level cfChIP enrichment returned by in_house_frag_counts()
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
GSEA_barplot <- function(gsea_res, #GSEA results returned by fgsea::fgsea()
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
basic_PCA <- function(counts, #Gene level cfChIP enrichment returned by in_house_frag_counts()
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
TSS_counts_selected_genes <- function(counts, #Gene level cfChIP enrichment returned by in_house_frag_counts()
                                      meta, #cfChIP metadata for each sample
                                      genes, #Gene SYMBOLS which should be plotted.
                                      samples, #List of samples, where the names of the 
                                      #list entries refer to the sample groups 
                                      #determining the x-value in the boxplot
                                      er = NULL, #Enrichment ratio cutoff
                                      lib_size = NULL, #Number of cfChIP fragments cutoff
                                      color = NULL, #Definition of colors to use,
                                      #otherwise viridis(begin = 0.3, end = 0.9) is used
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
#Plotting coverage around TSS for selected genes
coverage_selected_genes <- function(cov_arr, #The coverage relative to TSS for 
                                    #all genes and patients returned by gene_coverage_array()
                                    samples, #List of samples, where the names of the 
                                    #list entries refer to the sample groups 
                                    #determining the color in the plot.
                                    #If it is a nested list, the first layer
                                    #will determine the linetype
                                    meta, #cfChIP metadata for each sample
                                    genes, #Gene SYMBOLS which should be plotted.
                                    er = NULL, #Enrichment ratio cutoff
                                    lib_size = NULL, #Number of cfChIP fragments cutoff
                                    color = NULL, #Definition of colors to use, 
                                    #otherwise viridis(begin = 0.3, end = 0.9) is used
                                    concat = F, #If multiple genes are listed,
                                    #concat = T will plot the average enrichment,
                                    #if concat = F, individual plots will be created for each gene 
                                    title = NULL #Add plot title
                                    ){
    #Filtering samples
    all_samples <- unlist(samples)
    meta <- meta %>% 
        filter(sample %in% all_samples)
    meta <- filter_samples(meta = meta,
                           er = er,
                           lib_size = lib_size)
    #Evaluating whether coverage data for the genes exist
    mis_genes <- genes[!genes %in% dimnames(cov_arr)$gene]
    if(length(mis_genes) > 0){
        gene_list <- paste0(mis_genes, collapse = "," )
        message(paste0("The following gene is not included:\n",
                       gene_list))
        genes <- genes[genes %in% dimnames(cov_arr)$gene]
    }
    #Selecting samples
    sele_samples <- all_samples[all_samples %in% meta$sample]
    cov_sele <- cov_arr[,genes,sele_samples]
    list_of_lists <- is.list(samples[[1]])
    if(list_of_lists){
        #Generating plot if samples is a nested list
        sample_groups1 <- stack(samples[[1]])
        sample_groups1$timepoint <- names(samples)[1]
        sample_groups2 <- stack(samples[[2]])
        sample_groups2$timepoint <- names(samples)[2]
        sample_groups <- rbind(sample_groups1,sample_groups2)
        colnames(sample_groups)[1:2] <- c("individual", "group")
        if(length(genes) > 1 & !concat){
            #formatting data for multiple genes plotted individually
            df <- reshape2::melt(cov_sele, value.name = "coverage") %>% 
                left_join(sample_groups, by = "individual") %>% 
                mutate(gene = as.character(gene)) %>% 
                rename_AC(char = "gene") %>% 
                group_by(timepoint,group,gene,position) %>% 
                summarise(m_cov = mean(coverage)) %>% 
                ungroup() %>% 
                filter(position >= -2100 & position <= 2100) %>% 
                mutate(gene = factor(gene,levels = genes))
            tit <- " "
        } else{
            #formatting data for one gene or gene averages
            df <- reshape2::melt(cov_sele, value.name = "coverage") %>% 
                left_join(sample_groups, by = "individual") %>% 
                group_by(timepoint,group,position) %>% 
                summarise(m_cov = mean(coverage)) %>% 
                ungroup() %>% 
                filter(position >= -2160 & position <= 2160) 
            if(length(genes) == 1){
                df <- df %>% 
                    mutate(SYMBOL = genes) %>% 
                    rename_AC(char = "SYMBOL")
                tit <- unique(df$SYMBOL)
            }
        }
        #Defining colors
        if(is.null(color)){
            vals <- viridis(begin = 0.3,
                            end = 0.9,
                            n = length(samples[[1]]))
            names(vals) <- names(samples[[1]])
        }else{
            vals <- color
            names(vals) <- names(samples[[1]])
        }
    }
    else{
        #Formatting data if samples is not a nested list
        sample_groups <- stack(samples)
        colnames(sample_groups) <- c("individual", "group")
        if(length(genes) > 1 & !concat){
            #formatting data for multiple genes plotted individually
            df <- reshape2::melt(cov_sele, value.name = "coverage") %>% 
                left_join(sample_groups, by = "individual") %>% 
                mutate(gene = as.character(gene)) %>% 
                rename_AC(char = "gene") %>% 
                group_by(group,gene,position) %>% 
                summarise(m_cov = mean(coverage)) %>% 
                ungroup() %>% 
                filter(position >= -2100 & position <= 2100) %>% 
                mutate(gene = factor(gene,levels = genes))
            tit <- " "
        } else{
            #formatting data for one gene or gene averages
            df <- reshape2::melt(cov_sele, value.name = "coverage") %>% 
                left_join(sample_groups, by = "individual") %>% 
                group_by(group,position) %>% 
                summarise(m_cov = mean(coverage)) %>% 
                ungroup() %>% 
                filter(position >= -2160 & position <= 2160) 
            if(length(genes) == 1){
                df <- df %>% 
                    mutate(SYMBOL = genes) %>% 
                    rename_AC(char = "SYMBOL")
                tit <- unique(df$SYMBOL)
            }
        }
        #Defining colors
        if(is.null(color)){
            vals <- viridis(begin = 0.3,
                            end = 0.9,
                            n = length(samples))
            names(vals) <- names(samples)
        }else{
            vals <- color
            names(vals) <- names(samples)
        }
    }
    #Defining title
    if(length(genes) > 1 & is.null(title)){
        tit <- "Gene list"
    }
    if(!is.null(title)){
        tit <- title
    }
    if(list_of_lists){
        #Creating plot for samples as a nested list
        gg <- ggplot(df, 
                     aes(x = position, 
                         y = m_cov))+
            geom_smooth(aes(color = group,
                            linetype = timepoint,
                            group = interaction(group,timepoint)),
                        se = FALSE,
                        span = 0.2,
                        method = "loess",
                        formula = "y ~ x")+
            scale_linetype_discrete("Timepoint")+
            theme_classic()+
            guides(color = guide_legend(order = 1),
                   linetype = guide_legend(order = 2,
                                           override.aes = list(color = "black")))+
            theme(legend.position = c(0.85,0.70),
                  legend.background = element_rect(fill = alpha("white",0.8),
                                                   color = NA),
                  legend.key = element_rect(fill = NA))
    }else{
        #Creating plot if samples is not a nested list
        gg <- ggplot(df, 
                     aes(x = position, 
                         y = m_cov,
                         color = group))+
            theme_classic()+
            geom_smooth(se = F,
                        span = 0.2,
                        method = "loess",
                        formula = "y ~ x")+
            theme(legend.position = "right")
    }
    #Defining limits
    gg <- gg +
        scale_color_manual("Sample   ",values = vals)+
        scale_x_continuous(limits = c(-2160,2160),
                           breaks = c(-2000,0,2000),
                           labels = c("-2 kb", "TSS", "2 kb"),
                           expand = c(0,0))+
        scale_y_continuous(expand = c(0,0))+
        th
    #Adding labels
    if(length(genes) == 1 | concat){
        gg <- gg+
            labs(title = tit,
                 y = "H3K4me3 Coverage",
                 x = "")
    }
    else{
        gg <- gg+
            labs(y = paste0("H3K4me3 coverage"),
                 x = "")
    }
    if(!concat){
        #Creating plot for multiple individual genes
        gg <- gg +
            facet_wrap(~gene,
                       scales = "free_y",
                       nrow = 1)+
            scale_x_continuous(limits = c(-2160,2160),
                               breaks = c(-2000,0,2000),
                               labels = c("-2 kb", "TSS", "2 kb"),
                               expand = c(0,0))+
            theme(strip.background = element_blank(),
                  strip.text = element_text(face = "bold",
                                            size = 12),
                  legend.position = "top",
                  panel.background = element_rect(fill = "white",
                                                  color = "white"),
                  plot.background = element_rect(fill = "white",
                                                 color = "white"),
                  
                  panel.spacing = unit(25,"pt"))
    }
    return(gg)
}

# Surivival plots
survival_plot <- function(data, #Clinical data w. "divider" variable to separate patients into groups
                          survival, #OS or PFS to define which should be plotted
                          plot_title, #Add title
                          legend_title, #Add legend title
                          groups, #Names of groups in the legend
                          color, #Color for curves
                          p_x = 25 #p-value x-axis coordinate
                          ){
    #Calculating surivival fit for OS or PFS
    if(survival == "OS"){
        calc_fit <- survfit(Surv(OS,event_OS) ~divider,
                            data = data)
    }
    if(survival == "PFS"){
        calc_fit <- survfit(Surv(PFS,event_PFS) ~divider,
                            data = data)
    }
    #Calculating log-rank p-value
    surv_diff <- survdiff(formula = formula(calc_fit),
                          data = data)
    p_surv <- 1-pchisq(surv_diff$chisq,
                       df = length(surv_diff$n)-1)
    #Formatting p-value
    ps <- format_p(p_surv)
    #Creating plot
    gg <- ggsurvplot(data = data,
                     fit = calc_fit,
                     xlab = "Months",
                     ylab = "Probability of survival",
                     font.x = c(12,"bold","black"),
                     font.y = c(12,"bold","black"),
                     title = plot_title,
                     font.main = c(14,"bold","black"),
                     palette = color,
                     legend = c(0.8,0.8),
                     legend.title = legend_title,
                     legend.labs = groups,
                     font.legend = c(12,"plain","black"),
                     surv.median.line = "hv",
                     risk.table = TRUE,
                     risk.table.pos = "out",
                     tables.y.text = FALSE,
                     tables.title = "",
                     tables.theme = theme_cleantable(),
                     tables.height = 0.2,
                     censor.shape = "I",
                     font.tickslab = c(10))
    gg$plot <- gg$plot +
        annotate("text",
                 x = p_x,
                 label = paste0("p == ",ps),
                 parse = TRUE,
                 y = 0.62,
                 size = 4)
    return(gg)
}

# Plotting surivival based on clusters
survival_clusters <- function(clinical, #Clinical data
                              clusters, #List (length  = 2) of samples within each cluster
                              l_title, #Add legend title
                              color = NULL, #Define curve colors
                              p_x = 25, #p-value x-axis coordinate
                              get_coph = FALSE #If FALSE the plot will be generated, 
                              #if TRUE instead the surivival stats are calculated and exported
                              ){
    #Collecting clusters and samples
    cluster_df <- stack(clusters) %>% 
        dplyr::rename("patient" = "values") %>% 
        dplyr::rename("divider" = "ind") %>% 
        mutate(divider = divider)
    clinical_df <- clinical %>% 
        left_join(cluster_df, by = "patient")
    #Defining colors
    if(is.null(color)){
        vals <- viridis(begin = 0.3, 
                        end = 0.9,
                        n = length(clusters))
    }else{
        vals <- color
    }
    #Plotting OS of clusters
    gg_OS <- survival_plot(data = clinical_df,
                           survival = "OS",
                           plot_title = "OS",
                           legend_title = l_title,
                           groups = names(clusters),
                           color = vals,p_x = p_x)
    #Plotting PFS of clusters
    gg_PFS <- survival_plot(data = clinical_df,
                            survival = "PFS",
                            plot_title = "PFS",
                            legend_title = l_title,
                            groups = names(clusters),
                            color = vals,p_x = p_x)
    gg <- arrange_ggsurvplots(x = list(gg_PFS,gg_OS),
                              ncol = 2,
                              print = FALSE,
                              nrow = 1)
    #Calculating HR and stats on PFS and OS
    if(get_coph){
        meta_sele_BL <- meta %>% 
            filter(study_name == "Lung cancer") %>% 
            filter(sample_type == "Baseline") %>% 
            dplyr::select(patient,ctDNA) %>% 
            unique()
        clinical_df <- clinical_df %>% 
            left_join(meta_sele_BL, by = "patient")
        cox_pfs <- summary(coxph(Surv(PFS,event_PFS)
                                 ~ divider + patient_group, 
                                 data = clinical_df))
        cox_os <- summary(coxph(Surv(OS,event_OS) 
                                ~ divider + patient_group, 
                                data = clinical_df))
        fit_pfs <- summary(survfit(Surv(PFS, event_PFS)
                                   ~ divider,
                                   data = clinical_df))$table
        fit_os <- summary(survfit(Surv(OS, event_OS)
                                  ~ divider,
                                  data = clinical_df))$table
        return(list("HR PFS" = cox_pfs,
                    "HR OS" = cox_os,
                    "med PFS" = fit_pfs,
                    "med OS" = fit_os))
    }else{
        return(gg)
    }
}

#Survival based on gene enrichment
survival_gene_express <- function(counts, #Gene level cfChIP enrichment returned by in_house_frag_counts()
                                  meta, #cfChIP metadata for each sample
                                  SYMBOL, #Gene SYMBOL where the enrichment will separate patients
                                  samples, #Samples to be included in analysis
                                  clinical, #Clinical data
                                  er = NULL, #Enrichment ratio cutoff
                                  lib_size = NULL, #Number of cfChIP fragments cutoff
                                  expression = "continuous", #If "continuous" patients will 
                                  #be separated based on group median,
                                  #If "dynamic" gene dynamics will separate patients based on
                                  #ΔH3K4me3 enrichment > 0 or < 0
                                  color = NULL, #Define curve colors
                                  p_x = 25, #p-value x-axis coordinate
                                  get_coph = FALSE #If FALSE the plot will be generated, 
                                  #if TRUE instead the surivival stats are calculated and exported
                                  ){
    if(expression == "continuous"){
        #Dividing patients based on median enrichment
        divider_df <- dichotomize_patients_enrichment(counts = counts,
                                                      meta = meta,
                                                      samples = samples,
                                                      SYMBOL = SYMBOL,
                                                      er = er,
                                                      lib_size = lib_size) %>% 
            mutate(SYMBOL = SYMBOL) %>% 
            rename_AC(char = "SYMBOL")
        groups <- c("Low", "High")
        l_title <- paste0(unique(divider_df$SYMBOL), 
                          " enrichment")
    }
    if(expression == "dynamic"){
        #Dividing patients based on enrichment dynamics from BL to Tx
        divider_df <- dichotomize_patients_dynamics(counts = counts,
                                                    meta = meta,
                                                    SYMBOL = SYMBOL,
                                                    er = er,
                                                    lib_size = lib_size) %>% 
            mutate(SYMBOL = SYMBOL) %>% 
            rename_AC(char = "SYMBOL")
        groups <- c("Decrease", "Increase")
        l_title <- paste0(unique(divider_df$SYMBOL),
                          " dynamics")
    }
    #Defining colors
    if(is.null(color)){
        vals <- viridis(begin = 0.3,
                        end = 0.9,
                        n = length(groups))
        names(vals) <- names(groups)
    }else{
        vals <- color
        names(vals) <- names(groups)
    }
    #Creating plots
    df_surv <- divider_df %>% 
        left_join(clinical, by = "patient")
    gg_OS <- survival_plot(data = df_surv,
                           survival = "OS",
                           plot_title = "OS",
                           legend_title = l_title,
                           groups = groups,
                           color = vals,
                           p_x = p_x)
    gg_PFS <- survival_plot(data = df_surv,
                            survival = "PFS",
                            plot_title = "PFS",
                            legend_title = l_title,
                            groups = groups,
                            color = vals,
                            p_x = p_x)
    gg <- arrange_ggsurvplots(x = list(gg_PFS,gg_OS),
                              ncol = 2,
                              print = FALSE,
                              nrow = 1)
    #Calculating HR and stats on PFS and OS
    if(get_coph){
        meta_sele_BL <- meta %>% 
            filter(study_name == "Lung cancer") %>% 
            filter(sample_type == "Baseline") %>% 
            dplyr::select(patient,ctDNA) %>% 
            unique()
        df_surv <- df_surv %>% 
            left_join(meta_sele_BL, by = "patient")
        cox_pfs <- summary(coxph(Surv(PFS,event_PFS)
                                 ~ divider + patient_group, 
                                 data = df_surv))
        cox_os <- summary(coxph(Surv(OS,event_OS) 
                                ~ divider + patient_group, 
                                data = df_surv))
        fit_pfs <- summary(survfit(Surv(PFS, event_PFS)
                                   ~ divider,
                                   data = df_surv))$table
        fit_os <- summary(survfit(Surv(OS, event_OS)
                                  ~ divider,
                                  data = df_surv))$table
        return(list("HR PFS" = cox_pfs,
                    "HR OS" = cox_os,
                    "med PFS" = fit_pfs,
                    "med OS" = fit_os))
    }else{
        return(gg)
    }
}

# Heatmap of individual gene enrichments
heatmap_expression <- function(counts, #Gene level cfChIP enrichment returned by in_house_frag_counts()
                               meta, #cfChIP metadata for each sample
                               genes, #Gene SYMBOLS which should be plotted.
                               samples, #List of samples, where the names of the 
                               #list entries refer to the sample groups 
                               #determining how heatmap is separated
                               er = NULL, #Enrichment ratio cutoff
                               lib_size = NULL, #Number of cfChIP fragments cutoff
                               color = NULL, #Definition of colors to use,
                               #otherwise viridis(begin = 0.3, end = 0.9) is used
                               title = NULL #Add title to x-axis
                               ){
    #Filtering samples
    samples_df <- stack(samples) %>% 
        dplyr::rename("sample" = "values") %>% 
        dplyr::rename("group" = "ind")
    meta <- meta %>% 
        filter_samples(er = er,
                       lib_size = lib_size) %>% 
        filter(sample %in% samples_df$sample)
    #Selecting samples and genes
    counts_sele <- counts[genes,meta$sample] %>% 
        as.data.frame()
    counts_sele$SYMBOL <- rownames(counts_sele)
    rownames(counts_sele) <- NULL
    #Formatting enrichment data
    df <- counts_sele %>% 
        pivot_longer(-SYMBOL) %>% 
        dplyr::rename("sample" = "name") %>% 
        dplyr::rename("signal" = "value") %>% 
        left_join(samples_df, by = "sample")
    #Calculating z-scores
    df_z <- df %>% 
        group_by(SYMBOL) %>% 
        summarise(mean_signal = mean(signal),
                  sd_signal = sd(signal),
                  .groups = "drop_last")
    df <- df %>% 
        left_join(df_z, by = "SYMBOL") %>% 
        mutate(z_score = (signal-mean_signal)/sd_signal)
    #Defining plot orders
    df_sample_orders <- df %>% 
        group_by(group, sample) %>% 
        summarise(mean_z = mean(z_score),.groups = "drop_last") %>% 
        arrange(mean_z) %>% 
        ungroup()
    df_SYMBOL_orders <- df %>% 
        group_by(group,SYMBOL) %>% 
        summarise(mean_z = mean(z_score),.grpups = "drop_last") %>% 
        arrange(mean_z) %>% 
        ungroup()
    #Defining color
    if(is.null(color)){
        vals <- viridis(begin = 0.3, 
                        end = 0.9,
                        n = length(samples))
    }else{
        vals <- color
    }
    df <- df %>% 
        mutate(sample = factor(sample,
                               levels = df_sample_orders$sample)) %>% 
        mutate(SYMBOL = factor(SYMBOL, 
                               levels = rev(unique(df_SYMBOL_orders$SYMBOL))))
    #Creating plot
    gg <- ggplot(df, 
                 aes(x = SYMBOL,
                     y = sample))+
        geom_tile(aes(fill = z_score))+
        facet_grid(group ~ .,
                   scales = "free_y",
                   space = "free")+
        labs(x = title,
             y = "Patient")+
        scale_fill_gradient2("Z-score",
                             low = vals[1],
                             mid = "black",
                             high = vals[2],
                             midpoint = 0)+
        theme_void()+
        th+
        theme(axis.text.y = element_blank(),
              axis.text.x = element_text(size = 8,
                                         angle = 45,
                                         vjust = 1,
                                         hjust = 1),
              legend.position = "top",
              axis.line.y = element_blank(),
              axis.title.x = element_text(angle = 0),
              axis.title.y = element_text(angle = 90,
                                          vjust = 0.5),
              panel.border = element_rect(color = "black",
                                          fill = NA),
              panel.spacing = unit(5,"pt"),
              strip.background = element_rect(fill = "white"),
              panel.grid.major.x = element_blank(),
              strip.text = element_text(color = "black",
                                        face = "bold",
                                        size = 10,
                                        vjust = 1,
                                        angle = -90),
              plot.margin = margin(t = 0.5,
                                   r = 0.5,
                                   b = 0.5,
                                   l = 0.5,
                                   unit = "cm"))
    return(gg)
}

#Heatmap of gene enrichment dynamics
heatmap_dynamics <- function(diff_df, #Gene enrichment dynamics returned by calc_diff()
                             meta, #cfChIP metadata for each sample
                             genes, #Gene SYMBOLS which should be plotted.
                             samples, #List of samples, where the names of the 
                             #list entries refer to the sample groups 
                             #determining how heatmap is separated
                             er = NULL, #Enrichment ratio cutoff
                             lib_size = NULL, #Number of cfChIP fragments cutoff
                             color = NULL, #Definition of colors to use,
                             #otherwise viridis(begin = 0.3, end = 0.9) is used
                             title = NULL #Add title to x-axis
                             ){
    #Filtering samples
    meta <- meta %>% 
        filter(study_name == "Lung cancer") %>% 
        filter_samples(er = er,
                       lib_size = lib_size)
    sample_df <- stack(samples) %>% 
        dplyr::rename("patient" = "values") %>% 
        dplyr::rename("group" = "ind")
    #Selecting relevant genes and samples
    df <- diff_df %>% 
        dplyr::select(all_of(c("SYMBOL",meta$patient))) %>% 
        filter(SYMBOL %in% genes) %>% 
        rename_AC(char = "SYMBOL") %>% 
        pivot_longer(-SYMBOL) %>%
        dplyr::rename("patient" = "name") %>% 
        dplyr::rename("diff" = "value") %>% 
        left_join(sample_df, by ="patient")
    #Defining plot orders
    df_sample_orders <- df %>% 
        group_by(group,patient) %>% 
        summarise(mean_diff = mean(diff),.groups = "drop_last") %>% 
        arrange(mean_diff) %>% 
        ungroup()
    df_gene_orders <- df %>% 
        group_by(group,SYMBOL) %>% 
        summarise(mean_diff = mean(diff),.groups = "drop_last") %>% 
        arrange(mean_diff) %>% 
        ungroup()
    df <- df %>% 
        mutate(patient = factor(patient, 
                                levels = df_sample_orders$patient)) %>% 
        mutate(SYMBOL = factor(SYMBOL, 
                               levels = rev(unique(df_gene_orders$SYMBOL))))
    #Defining colors
    if(is.null(color)){
        vals <- viridis(begin = 0.3, 
                        end = 0.9,
                        n = length(samples))
    }else{
        vals <- color
    }
    #Creating plot
    gg <- ggplot(df, 
                 aes(y = patient,
                     x = SYMBOL))+
        geom_tile(aes(fill = diff))+
        facet_grid(group ~ .,
                   scales = "free_y",
                   space = "free")+
        labs(title = "",
             x = title,
             y = "Patient")+
        scale_fill_gradient2(expression(bold(Delta * "H3K4me3 enrichment        ")),
                             low = vals[1],
                             mid = "black",
                             high = vals[2],
                             midpoint = 0,
                             breaks = c(-1,0,1),
                             guide = guide_colorbar(barwidth = unit(50,
                                                                    "pt")))+
        theme_void()+
        th+
        theme(axis.text.y = element_blank(),
              axis.text.x = element_text(angle = 45,
                                         vjust= 1, 
                                         hjust = 1),
              legend.position = "top",
              axis.line.y = element_blank(),
              axis.title.x = element_text(angle = 0),
              axis.title.y = element_text(angle = 90,
                                          vjust= 0.5),
              panel.border = element_rect(color = "black",
                                          fill = NA),
              panel.spacing = unit(5,"pt"),
              strip.background = element_rect(fill = "white"),
              panel.grid.major.x = element_blank(),
              strip.text = element_text(color = "black",
                                        face = "bold",
                                        size = 10,
                                        vjust = 1,
                                        angle = -90))
    return(gg)
}

#Plotting differences in gene enrichment
gene_dynamic_plot <- function(diff_df, #Gene enrichment dynamics returned by calc_diff()
                              meta, #cfChIP metadata for each sample
                              genes, #Gene SYMBOLS which should be plotted.
                              er = NULL, #Enrichment ratio cutoff
                              lib_size = NULL, #Number of cfChIP fragments cutoff
                              color = NULL, #Definition of colors to use,
                              #otherwise viridis(begin = 0.3, end = 0.9) is used
                              title = NULL, #Add title to x-axis
                              patient_groups = NULL, #If NULL, then mR is compared to mNR
                              #Otherwise patient_groups = <data.frame> can determine which
                              #groups are compared
                              tit = NULL,
                              concat = TRUE){
    #Filtering samples
    meta <- meta %>% 
        filter(study_name == "Lung cancer") %>% 
        filter_samples(er = er,
                       lib_size = lib_size) %>% 
        dplyr::select(patient,patient_group)
    c_meta <- meta %>% 
        group_by(patient) %>% 
        summarise(n = n())
    meta <- meta %>% 
        left_join(c_meta, by = "patient") %>% 
        filter(n == 2) %>% 
        unique()
    if(length(genes) == 1){
        #Formatting data for a single gene
        df <- diff_df %>% 
            dplyr::select(all_of(c("SYMBOL",meta$patient))) %>% 
            filter(SYMBOL %in% genes) %>% 
            rename_AC(char = "SYMBOL") %>% 
            pivot_longer(-SYMBOL) %>% 
            dplyr::rename("patient" = "name") %>% 
            dplyr::rename("diff_signal" = "value") %>% 
            left_join(meta, by = "patient") %>% 
            unique() %>% 
            mutate(patient_group = ifelse(patient_group == "responder",
                                          "mR",
                                          "mNR"))
        tit <- unique(df$SYMBOL)
        ylab <- expression(bold(Delta * "H3K4me3 enrichment"))
    }else{
        #Formatting data for a multiple genes
        df <- diff_df %>% 
            dplyr::select(all_of(c("SYMBOL",meta$patient))) %>% 
            filter(SYMBOL %in% genes)
        if(concat){
            #Calculating average across genes
            m_vals <- df %>% 
                dplyr::select(-SYMBOL) %>% 
                t %>% 
                rowMeans()
            df <- data.frame(patient = names(m_vals),
                             diff_signal = m_vals) 
            ylab <- expression(bold("Mean" ~ Delta * "H3K4me3 enrichment"))
        }else{
            #Foratting data to plot individual genes
            df <- df %>% 
                pivot_longer(-SYMBOL) %>% 
                dplyr::rename("patient" = "name") %>% 
                dplyr::rename("diff_signal" = "value")
            ylab <- expression(bold(Delta * "H3K4me3 enrichment"))
        }
        df <- df %>% 
            left_join(meta, by = "patient") %>% 
            unique() %>% 
            mutate(patient_group = ifelse(patient_group == "responder",
                                          "mR",
                                          "mNR"))
    }
    if(!is.null(patient_groups)){
        #If patient_groups info is added the info is incorporated
        df <- df %>% 
            dplyr::select(-patient_group) %>% 
            left_join(patient_groups, by = "patient")
    }
    #Defining colors
    if(is.null(color)){
        vals <- viridis(begin = 0.2, end = 0.8,n = 2)
    }else{
        vals <- color
    }
    #Defining comparisons
    compares <- combn(as.character(unique(df$patient_group)),
                      m = 2,
                      simplify = F)
    #Performing statistics
    if(concat){
        stat_test <- df %>% 
            t_test(diff_signal ~ patient_group,
                   comparisons = compares,
                   paired = FALSE) %>% 
            add_xy_position(x = "patient_group",
                            data = df,
                            formula = diff_signal ~ patient_group,
                            scales = "free_y",
                            step.increase = 0.12) 
    }else{
        stat_test <- df %>% 
            group_by(SYMBOL) %>% 
            t_test(diff_signal ~ patient_group,
                   comparisons = compares,
                   paired = FALSE) %>% 
            add_xy_position(x = "patient_group",
                            data = df,
                            formula = diff_signal ~ patient_group,
                            scales = "free_y",
                            step.increase = 0.12)
    }
    #Creating plot
    gg <- ggplot(df,
                 aes(x = patient_group,
                     y = diff_signal,
                     fill = patient_group))+
        geom_boxplot(outlier.alpha = 0)+
        theme_classic()+
        g_brac(stat_test)+
        scale_y_continuous(expand = c(0.2,0))+
        th+
        theme(legend.position = "none")+
        geom_hline(yintercept = 0,
                   linetype = "dashed",
                   linewidth = 1)+
        scale_fill_manual(values = vals)
    if(concat | length(genes) == 1){
        gg <- gg+
            labs(title = tit,
                 y = ylab,
                 x = "")
    }else{
        gg <- gg+
            facet_wrap(~SYMBOL,
                       nrow = 1, 
                       scales = "free")+
            labs(y = ylab,
                 x = "")+
            scale_y_continuous(expand = c(0.2,0))+
            theme(strip.background = element_blank(),
                  strip.text = element_text(face = "bold",
                                            size = 12),
                  panel.background = element_rect(fill = "white",
                                                  color = "white"),
                  plot.background = element_rect(fill = "white",
                                                 color = "white"))
    }
    return(gg)
}

#Plotting gene expression percentiles
percentile_violin <- function(gte, #GTEx gene expression percentiles returned by calc_individual_percentile()
                              ruv, #RUV-III gene expression percentiles returned by calc_individual_percentile()
                              genes, #Gene SYMBOLS which should be plotted.
                              color = NULL #Definition of colors to use,
                              #otherwise viridis(begin = 0.3, end = 0.9) is used
                              ){
    #Collecting data
    gte_data <- gte %>% 
        filter(SYMBOL == genes) %>% 
        pivot_longer(-SYMBOL) %>% 
        mutate(group = "Whole blood")
    ruv_data <- ruv %>% 
        filter(SYMBOL == genes) %>% 
        pivot_longer(-SYMBOL) %>% 
        mutate(group = "NSCLC")
    df <- rbind(gte_data,ruv_data) 
    #Defining color
    if(is.null(color)){
        vals <- viridis(begin = 0.3, 
                        end = 0.9,
                        n = length(unique(df$group)))
    }else{
        vals <- color
    }
    #Establish comparisons
    names(vals) <- c("NSCLC", "Whole blood")
    compares <- combn(unique(df$group),m = 2,
                      simplify = F)
    stat_test <- df %>% 
        wilcox_test(value ~ group,
                    comparisons = compares,
                    paired = FALSE) %>% 
        add_xy_position(x = "group",
                        data = df,
                        formula = value ~ group,
                        scales = "free_y",
                        step.increase = 0.12)
    
    df <- df %>% 
        mutate(group = factor(group, levels = c("Whole blood", 
                                                "NSCLC")))
    #Creating plot
    gg <- ggplot(df,aes(x = group,
                        y = value,
                        fill = group))+
        geom_violin(scale = "width")+
        scale_fill_manual("Tissue",values = vals)+
        geom_jitter(width = 0.2,
                    alpha = 0.2)+
        g_brac(stat_test)+
        labs(title = genes,
             y = "Rank",
             x = "")+
        stat_compare_means(comparisons = compares,
                           method = "wilcox.test")+
        scale_y_continuous(limits = c(0,16000))+
        theme_classic()+
        th+
        theme(legend.position = "none")
    return(gg)
}

#Plotting paired enrichment for BL and Tx
interaction_plot <- function(counts, #Gene level cfChIP enrichment returned by in_house_frag_counts()
                             meta, #cfChIP metadata for each sample
                             SYMBOL, #Gene SYMBOLS which should be plotted.
                             er = NULL, #Enrichment ratio cutoff
                             lib_size = NULL #Number of cfChIP fragments cutoff
                             ){
    #Filtering samples
    meta <- meta %>% 
        filter(study_name == "Lung cancer") %>% 
        filter_samples(er = er,
                       lib_size = lib_size)
    c_meta <- meta %>% 
        group_by(patient) %>% 
        summarise(n = n())
    meta <- meta %>% 
        left_join(c_meta, by = "patient") %>% 
        filter(n == 2)
    #Selecting genes and samples
    counts_sele <- counts[SYMBOL,meta$sample]
    if(length(SYMBOL) == 1){
        #Formatting data for a single gene
        df <- as.data.frame(counts_sele)
        colnames(df) <- "signal"
        df$sample <- rownames(df)
        rownames(df) <- NULL
        df <- df %>% 
            mutate(SYMBOL = SYMBOL) %>% 
            rename_AC(char = "SYMBOL")
        char <- unique(df$SYMBOL)
    }else{
        #Formatting data for multiple genes
        df <- counts_sele %>% 
            as.data.frame()
        df$SYMBOL <- rownames(df)
        df <- df %>% 
            pivot_longer(-SYMBOL) %>% 
            dplyr::rename("sample" = "name") %>% 
            dplyr::rename("signal" = "value")
    }
    meta_sele <- meta %>% 
        dplyr::select(sample, patient,sample_type,patient_group)
    df <- df %>% 
        left_join(meta_sele, by = "sample") %>% 
        mutate(patient_group = ifelse(patient_group == "responder",
                                      "mR",
                                      "mNR")) %>% 
        mutate(sample_type = ifelse(sample_type == "Baseline","BL","Tx"))
    #Creating plot
    gg <- ggplot(df,
                 aes(x = sample_type,
                     y = signal,
                     group = patient,
                     color = patient_group))+
        geom_line(alpha = 0.6)+
        geom_point(size = 2)+
        scale_color_viridis_d("Patient group    ", 
                              begin = 0.2, 
                              end = 0.8)+
        theme_classic()+
        th+
        theme(legend.position = "top")
    if(length(SYMBOL) == 1){
        gg <- gg+
            labs(title = char,
                 y = "H3K4me3 enrichment",
                 x = "Timepoint")+
            guides(color = guide_legend(ncol = 1, 
                                        byrow = TRUE))
    }
    else{
        gg <- gg+
            facet_wrap(~SYMBOL,
                       nrow = 1,
                       scales = "free")+
            labs(y = "H3K4me3 enrichment",
                 x = "Timepoint")+
            theme(strip.background = element_blank(),
                  strip.text = element_text(face = "bold",
                                            size = 12),
                  panel.background = element_rect(fill = "white",
                                                  color = "white"),
                  plot.background = element_rect(fill = "white",
                                                 color = "white"))
    }
    return(gg)
}

#Plotting gene annotations
peak_stats <- function(peaks, #Peak annotations across all samples
                       samples #Samples to include in plot
                       ){
    samples_df <- stack(samples) %>% 
        dplyr::rename("sample" = "values") %>% 
        dplyr::rename("tit" = "ind")
    #Filtering samples and calculating peak stats
    df <- peaks %>% 
        filter(sample %in% samples_df$sample) %>% 
        mutate(annotation = new_annotation(annotation)) %>% 
        group_by(sample,annotation) %>% 
        summarise(n = n()) %>% 
        ungroup() %>% 
        group_by(sample) %>% 
        reframe(total = sum(n),
                Feature = annotation,
                n = n,
                Frequency = n/total * 100) %>% 
        left_join(samples_df, by = "sample") %>% 
        arrange(desc(Frequency))
    df <- df %>% 
        mutate(tit = factor(tit, 
                            levels = samples_df$tit)) %>% 
        mutate(Feature = factor(Feature, 
                                levels = unique(df$Feature)))
    #Creating plot
    gg <- ggplot(df, 
                 aes(x = tit,
                     y = Frequency))+
        geom_bar(aes(fill = Feature),
                 stat = "identity")+
        scale_fill_viridis_d(name = "Annotation",
                             begin = 0.2,
                             end = 0.9, 
                             option = "F")+
        labs(x = "Sample",
             y = "% Of peaks",
             title = "H3K4me3 cfChIP peaks")+
        geom_text(aes(label = ifelse(Frequency > 3,
                                     paste0(n),
                                     ""),
                      group = Feature),
                  size = 3,
                  color = "white",
                  position = position_stack(vjust = 0.5))+
        geom_text(aes(label = total,
                      x = tit,
                      y = 105),
                  size = 4,
                  vjust = 0.5,
                  color = "#751F58FF")+
        scale_y_continuous(limits = c(0,110),
                           breaks = seq(0,100,10),
                           expand = c(0,0))+
        theme_classic()+
        th+
        theme(legend.position = "right",
              legend.text = element_text(size = 8.5),
              axis.text.x = element_text(angle = 30,
                                         vjust= 1, 
                                         hjust = 1),
              legend.title = element_text(size = 10,
                                          face = "bold"))
    return(gg)
}