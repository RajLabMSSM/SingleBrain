suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readxl))


# For cell spc
df_coloc <- read_xlsx('~/all_COLOC_results_merged_H4_0.5_with_LD.xlsx')
df_dict <- df_coloc %>% dplyr::filter(PP.H4.abf>=0.80)

df_dict$feature <- sapply(strsplit(df_dict$feature, '\\.'), '[[', 1)
df_dict <- df_dict %>% drop_na()

df_sum <- read_tsv('~/sumstat_coloc_H4_0.8.tsv')
df_sum$labels <-paste0(df_sum$celltype,'\n','Beta: ',df_sum$beta,'\n','P: ',df_sum$pvalues)

snp_list <- df_dict$SNP %>% unique()

for(snp in snp_list){
  df_dict_tmp <-df_dict %>% dplyr::filter(SNP==snp) %>% distinct(SNP, Symbol, feature, ID)
  gene_list <- df_dict_tmp$Symbol %>% unique()

  for(genename in gene_list){

    df_dict_tmp2 <- df_dict_tmp %>% dplyr::filter(Symbol==genename) 
    geneid <- df_dict_tmp2$feature
    message(paste0(snp,'-',genename,'-',geneid))
    
    #Gene expression
    df_Ast_tmp <- df_Ast %>% dplyr::filter(Gene==geneid) %>% dplyr::select(-Gene)
    df_End_tmp <- df_End %>% dplyr::filter(Gene==geneid) %>% dplyr::select(-Gene)
    df_Ext_tmp <- df_Ext %>% dplyr::filter(Gene==geneid) %>% dplyr::select(-Gene)
    df_IN_tmp <- df_IN %>% dplyr::filter(Gene==geneid) %>% dplyr::select(-Gene)
    df_MG_tmp <- df_MG %>% dplyr::filter(Gene==geneid) %>% dplyr::select(-Gene)
    df_OD_tmp <- df_OD %>% dplyr::filter(Gene==geneid) %>% dplyr::select(-Gene)
    df_OPC_tmp <- df_OPC %>% dplyr::filter(Gene==geneid) %>% dplyr::select(-Gene)
    
    df_Ast_tmp = setNames(data.frame(t(df_Ast_tmp[,-1])), 'exp') #df_Ast_tmp[,1])
    df_End_tmp = setNames(data.frame(t(df_End_tmp[,-1])), 'exp') #df_End_tmp[,1])
    df_Ext_tmp = setNames(data.frame(t(df_Ext_tmp[,-1])), 'exp') #df_Ext_tmp[,1])
    df_IN_tmp = setNames(data.frame(t(df_IN_tmp[,-1])), 'exp') #df_IN_tmp[,1])
    df_MG_tmp = setNames(data.frame(t(df_MG_tmp[,-1])), 'exp') #df_MG_tmp[,1])
    df_OD_tmp = setNames(data.frame(t(df_OD_tmp[,-1])), 'exp') #df_OD_tmp[,1])
    df_OPC_tmp = setNames(data.frame(t(df_OPC_tmp[,-1])), 'exp') #df_OPC_tmp[,1])
    
    df_Ast_tmp$sample_id <- rownames(df_Ast_tmp)
    rownames(df_Ast_tmp) <- NULL
    df_End_tmp$sample_id <- rownames(df_End_tmp)
    rownames(df_End_tmp) <- NULL
    df_Ext_tmp$sample_id <- rownames(df_Ext_tmp)
    rownames(df_Ext_tmp) <- NULL
    df_IN_tmp$sample_id <- rownames(df_IN_tmp)
    rownames(df_IN_tmp) <- NULL
    df_MG_tmp$sample_id <- rownames(df_MG_tmp)
    rownames(df_MG_tmp) <- NULL
    df_OD_tmp$sample_id <- rownames(df_OD_tmp)
    rownames(df_OD_tmp) <- NULL
    df_OPC_tmp$sample_id <- rownames(df_OPC_tmp)
    rownames(df_OPC_tmp) <- NULL
    
    df_plot <- data.frame()
    for(celltype in cell_list){
      
      df_sum_tmp <- df_sum[df_sum$celltype==celltype,]
      df_sum_tmp <- df_sum_tmp[df_sum_tmp$QTL_Ensembl==geneid,]
      df_sum_tmp <- df_sum_tmp[df_sum_tmp$variant_id==snp,] # For GWAS
      
      if(nrow(df_sum_tmp)==0){ next }
      ref<- df_sum_tmp$ref
      alt<- df_sum_tmp$alt
      
      refhom <- paste0(ref,'/',ref)
      het <- paste0(ref,'/',alt)
      althom <- paste0(alt,'/',alt)
      
      df_geno_tmp <- df_geno %>% dplyr::filter(ID==snp)
      df_geno_tmp = setNames(data.frame(t(df_geno_tmp[,-1])), 'geno')
      df_geno_tmp <- df_geno_tmp %>% drop_na(geno)
      df_geno_tmp$sample_id <- rownames(df_geno_tmp)
      rownames(df_geno_tmp) <- NULL
      
      df_geno_tmp$genoanno <- df_geno_tmp$geno
      df_geno_tmp$genoanno <- gsub(2, althom,df_geno_tmp$genoanno)
      df_geno_tmp$genoanno <- gsub(1, het,df_geno_tmp$genoanno)
      df_geno_tmp$genoanno <- gsub(0, refhom,df_geno_tmp$genoanno)

      if(nrow(df_geno_tmp)==0){ next }
      if(celltype=='Ast'){
        df_celltype_tmp <- merge(df_geno_tmp, df_Ast_tmp, all.x = TRUE)
      }else if(celltype=='End'){
        df_celltype_tmp <- merge(df_geno_tmp, df_End_tmp, all.x = TRUE)
      }else if(celltype=='Ext'){
        df_celltype_tmp <- merge(df_geno_tmp, df_Ext_tmp, all.x = TRUE)
      }else if(celltype=='IN'){
        df_celltype_tmp <- merge(df_geno_tmp, df_IN_tmp, all.x = TRUE)
      }else if(celltype=='OD'){
        df_celltype_tmp <- merge(df_geno_tmp, df_OD_tmp, all.x = TRUE)
      }else if(celltype=='OPC'){
        df_celltype_tmp <- merge(df_geno_tmp, df_OPC_tmp, all.x = TRUE)
      }else if(celltype=='MG'){
        df_celltype_tmp <- merge(df_geno_tmp, df_MG_tmp, all.x = TRUE)
      }
      
      df_celltype_tmp$celltype <- celltype
      
      df_plot <- rbind(df_plot, df_celltype_tmp )
      
    }
    
    df_plot_fin <- df_plot  %>% drop_na()
    
    df_sum_tmp2 <- df_sum %>% dplyr::filter(QTL_Ensembl==geneid) %>%  dplyr::filter(variant_id==snp) %>%  
      dplyr::select(celltype, labels)
    df_plot_fin <- merge(df_plot_fin, df_sum_tmp2)
    
    df_plot_fin$geno <- factor(df_plot_fin$geno , levels = c(0, 1, 2))
    df_plot_fin$genoanno <- factor(df_plot_fin$genoanno,
                                   levels =c(refhom, het, althom)
    )
    
    label_names <- setNames(as.list(df_sum_tmp2$labels), nm = df_sum_tmp2$celltype)
    label_names2 <- c(label_names['Ext'],
                      label_names['IN'],
                      label_names['Ast'],
                      label_names['OD'],
                      label_names['OPC'],
                      label_names['MG'],
                      label_names['End']
    )
    label_level <- label_names2[!sapply(label_names2,is.null)]
    df_plot_fin$labels <- factor(df_plot_fin$labels ,
                                 levels = label_level
    )
    
    if(nrow(df_plot_fin)==0){next}
    p_box <-   ggplot(df_plot_fin, aes_string( y='exp', x='genoanno',color='celltype',alpha=0)) +
      geom_jitter(shape=16, 
                  na.rm = TRUE,
                  width = 0.2, height = 0,
                  alpha = 0.5,
                  # position=position_jitter(0.2)
      ) +
      scale_color_manual(values = colorset) +
      scale_fill_manual(values = colorset) +
      facet_grid(. ~ labels,) +
      geom_boxplot(outlier.shape = NA, color='black') + 
      ylab(paste0(genename," Residual expression")) +
      ggtitle(paste0(genename,
                     " - ", snp 
      )) +
      theme_classic() +
      theme(title = element_text(colour = "black", size=12, face = 'bold'),
            axis.text.x = element_text(colour = "black", size=12),
            axis.text.y = element_text(colour = "black",size=12),
            axis.title.x = element_blank(),
            axis.title.y = element_text(colour = "black",size=12),
            panel.spacing = unit(0, "mm"),
            strip.placement = "outside",
            strip.text = element_text(colour = "black",  size=10),
            strip.background = element_blank(),
            legend.position = "none",
      ) 
    
    show(p_box)
    cell_counts <- length(unique(df_plot_fin$celltype))    
    ggsave(filename, plot = p_box, units=c('in'),
           width = 1.2*cell_counts, height = 3, dpi = 600,
           limitsize = FALSE
    )
  }
}




