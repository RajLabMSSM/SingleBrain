
# Load the package 
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(splitstackshape))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v79))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(wesanderson))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readxl))


# Load COLOC & MR
COLOC_total <- read_tsv(paste0(path_dir,"/all_COLOC_results_merged_H4_0.5_with_LD.tsv.gz"))

MR_total <- read_tsv(paste0(path_dir,"/all_MR_results_merged_H4_1_snp_wFDR.tsv.gz"))
MR_total <- MR_total %>% dplyr::filter(grepl("Inverse variance weighted", method)) 

## Extract specific GWAS
COLOC_result <- COLOC_total %>% dplyr::filter(GWAS==gwas_name)

AD_total_colocs <- COLOC_result

### Convert gene id to gene symbol
genes <- AD_total_colocs$QTL_Gene[grepl('ENSG0', AD_total_colocs$QTL_Gene, fixed=TRUE)]
genes <- sapply(strsplit(genes, '\\.'), '[[', 1)
AD_total_colocs$QTL_Ensembl[is.na(AD_total_colocs$QTL_Ensembl)] = AD_total_colocs$QTL_Gene[is.na(AD_total_colocs$QTL_Ensembl)]
AD_total_colocs$QTL_Ensembl <- sapply(strsplit(AD_total_colocs$QTL_Ensembl, '\\.'), '[[', 1)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
colnames(geneIDs1) = c('new_QTL_Gene', 'QTL_Ensembl')
AD_total_colocs <- left_join(AD_total_colocs, geneIDs1, by='QTL_Ensembl')
AD_total_colocs$new_QTL_Gene[is.na(AD_total_colocs$new_QTL_Gene)] = AD_total_colocs$QTL_Gene[is.na(AD_total_colocs$new_QTL_Gene)]
AD_total_colocs <- AD_total_colocs[, c("GWAS","locus", 'GWAS_P',"PP.H4.abf","new_QTL_Gene","QTL", 'celltype')]
AD_total_colocs$type <- AD_total_colocs$QTL

AD_total_colocs$category = '[0.5, 0.7]'
AD_total_colocs$category[AD_total_colocs$PP.H4.abf >= 0.7] = '[0.7, 0.9]'
AD_total_colocs$category[AD_total_colocs$PP.H4.abf >= 0.9] = '[0.9, 1]'

AD_total_colocs <- AD_total_colocs[which(AD_total_colocs$QTL %in% snrna_qtl ),]

AD_total_colocs$category <- factor(AD_total_colocs$category , levels = c('[0.5, 0.7]', '[0.7, 0.9]','[0.9, 1]'))

### Make a summary count of coloc
AD_total_colocs_summary <- AD_total_colocs %>% 
  arrange(QTL, desc(PP.H4.abf)) %>% 
  group_by(QTL, celltype) %>% 
  distinct(locus, .keep_all = TRUE) %>% count(category)

AD_total_by_qtl <- AD_total_colocs_summary %>% group_by(QTL, celltype) %>% summarise(total=sum(n))

AD_total_colocs_summary$QTL <- factor(AD_total_colocs_summary$QTL , levels = dataset_coloc)
AD_total_colocs_summary$celltype <- factor(AD_total_colocs_summary$celltype , levels = dataset_celltype)
AD_total_by_qtl$celltype <- factor(AD_total_by_qtl$celltype , levels = dataset_celltype)

AD_total_colocs$category <- factor(AD_total_colocs$category , levels = c('[0.5, 0.7]', '[0.7, 0.9]','[0.9, 1]'))

AD_total_colocs$type <- factor(AD_total_colocs$type , levels = dataset_coloc)
AD_total_colocs$QTL <- factor(AD_total_colocs$QTL , levels = dataset_coloc)

AD_total_colocs$PP.H4.abf <- as.numeric(AD_total_colocs$PP.H4.abf)
AD_total_colocs$PP.H4.abf <- AD_total_colocs$PP.H4.abf <- round(AD_total_colocs$PP.H4.abf, digits=2)
coloc_combined_melted <- melt(AD_total_colocs, ID='locus')
coloc_combined_melted_filt <-  coloc_combined_melted[order(coloc_combined_melted[,'new_QTL_Gene'],
                                                           coloc_combined_melted[,'type'], -coloc_combined_melted[,'value']),]
df <- distinct(coloc_combined_melted_filt, locus, new_QTL_Gene, type, .keep_all= TRUE)
df$celltype <- factor(df$celltype , levels = dataset_celltype)

### Extract eQTL PP4 > 0.8 
subsetted <- df
subsetted <- subsetted %>% dplyr::filter(value >= 0.8 )
subsetted$value <- 1

subsetted <- subsetted %>% arrange(QTL, desc(value)) %>% 
  group_by(QTL, locus) %>% 
  distinct(QTL, .keep_all = TRUE) 

### Load MR result and filter the significant result
MR_result <- MR_total %>% dplyr::filter(GWAS==gwas_name)
MR_result$FDR <- p.adjust(MR_result$pval, method = "BH" )
MR_result <- MR_result %>% dplyr::filter(FDR < 0.05)

MR_result <- MR_result[which(MR_result$QTL %in% snrna_qtl ),]
AD_total_mrs <- MR_result
AD_total_mrs <- AD_total_mrs[, c("GWAS","locus", 'GWAS_SNP','GWAS_P',"outcome",'exposure_Gene',
                                 "QTL", 'b','se','celltype')]
AD_total_mrs$type <- AD_total_mrs$QTL
AD_total_mrs$category = 'positive'
AD_total_mrs$category[AD_total_mrs$b < 0 ] ='negative'

names(AD_total_mrs)[names(AD_total_mrs) == 'b'] <- 'Beta'
names(AD_total_mrs)[names(AD_total_mrs) == 'se'] <- 'Se'
names(AD_total_mrs)[names(AD_total_mrs) == 'exposure_Gene'] <- 'new_QTL_Gene'

AD_total_mrs$Beta <- as.numeric(AD_total_mrs$Beta)
AD_total_mrs$Beta <- AD_total_mrs$Beta <- round(AD_total_mrs$Beta, digits=2)
AD_total_mrs$value <- AD_total_mrs$Beta
AD_total_mrs$Beta <- abs(AD_total_mrs$Beta)

### Format MR and COLOC to merge
subsetted2 <- AD_total_mrs %>% dplyr::select(GWAS, locus, new_QTL_Gene, QTL, celltype, value)
subsetted2$method <- 'MR_positive'
subsetted2[(subsetted2$value<0),'method'] = 'MR_negative'
subsetted2$value <- 1
subsetted2$type <- paste0('MR_', subsetted2$QTL)
subsetted2$id <- paste0(subsetted2$QTL,'_',subsetted2$locus,'_',subsetted2$new_QTL_Gene)

subsetted <- subsetted %>% dplyr::select(GWAS, locus, new_QTL_Gene, QTL, celltype, value)
subsetted$type <- paste0('COLOC_', subsetted$QTL)
subsetted$method <- 'COLOC'
subsetted$id <- paste0(subsetted$QTL,'_',subsetted$locus,'_',subsetted$new_QTL_Gene)
subsetted3 <- rbind(subsetted, subsetted2)

dataset_type <- c('COLOC_MiGA', 'MR_MiGA',
                  'COLOC_MG', 'MR_MG',
                  'COLOC_Ext','MR_Ext',
                  'COLOC_IN','MR_IN',
                  'COLOC_OD','MR_OD',
                  'COLOC_Ast','MR_Ast',
                  'COLOC_OPC','MR_OPC',
                  'COLOC_End','MR_End'
)

dataset_coloc <- c('MiGA',
                   'MG', 
                   'Ext',
                   'IN',
                   'OD',
                   'Ast',
                   'OPC',
                   'End'
)
dataset_celltype <- c('Microglia','Excitatory Neurons','Inhibitory Neurons',
                      'Oligodendrocytes','Astrocytes','OPCs',
                      'Endothelial cells')
subsetted3$celltype <- factor(subsetted3$celltype , levels = dataset_celltype)

subsetted4 <- rbind(subsetted2[which(subsetted2$id %in% subsetted$id ),],subsetted)
subsetted4$type <- factor(subsetted4$type , levels = dataset_type)
subsetted4$QTL <- factor(subsetted4$QTL , levels = dataset_coloc)

table(subsetted4$type)

p1 <- subsetted4 %>% ggplot( aes(x=type, y=new_QTL_Gene, color=celltype,  fill=celltype , shape=method)) + 
  geom_point(aes(size = value))+
  scale_shape_manual(values = c(16, 25, 24))+
  scale_fill_manual(values=colorset) +
  scale_colour_manual(values=colorset) +
  theme_classic() +
  facet_grid(locus ~ QTL, scales = "free", space = "free", switch = 'x' ) +
  ylab('') +
  xlab('') +
  theme(title = element_text( face='bold', colour = "black", size=10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "italic", colour = "black",size=10),
        axis.ticks.x = element_blank(),
        axis.ticks = element_line(colour = "black", size=0.25),
        axis.line = element_line(colour = "black", size=0.25), 
        strip.text.y = element_text(angle = 0, face='italic', colour = "black", hjust = 0, size=10),
        strip.text.x = element_text(colour = "black", size=10),
        strip.placement = "bottom",
        strip.background = element_blank(),
        legend.position = "none",
        legend.box.margin = margin(c(0,0,0,0)),
        legend.margin = margin(c(0,0,0,0)),
        panel.spacing.y = unit(x = 0, units = "points"), 
        panel.border = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
  ) +
  ggtitle(disease) 

show(p1)

ggsave(filename, plot = p1, 
       width=6.5 ,height=13, 
       dpi=600)


## COLOC bar plot from coloc summary
AD_total_colocs_summary$QTL <- factor(AD_total_colocs_summary$QTL , levels = dataset_coloc)
AD_total_colocs_summary$celltype <- factor(AD_total_colocs_summary$celltype , levels = dataset_celltype)
AD_total_by_qtl$celltype <- factor(AD_total_by_qtl$celltype , levels = dataset_celltype)

AD_total_colocs$category <- factor(AD_total_colocs$category , levels = c('[0.5, 0.7]', '[0.7, 0.9]','[0.9, 1]'))
AD_total_colocs$type <- factor(AD_total_colocs$type , levels = dataset_coloc)
AD_total_colocs$QTL <- factor(AD_total_colocs$QTL , levels = dataset_coloc)

AD_total_colocs_summary <- AD_total_colocs_summary %>% dplyr::filter(category=='[0.9, 1]')
AD_total_by_qtl <- AD_total_colocs_summary %>% group_by(QTL, celltype) %>% summarise(total=sum(n))
AD_total_by_qtl$celltype <- factor(AD_total_by_qtl$celltype , levels = dataset_celltype)

barplot <- ggplot(AD_total_colocs_summary, aes(y=n, x=QTL, fill=celltype)) + #alpha=category, 
  geom_bar(position="stack", stat="identity") + facet_grid(. ~celltype, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = colorset) + 
  scale_alpha_discrete(name = "Category", range = c(0.3, 1)) + 
  theme_classic() + ylab('# GWAS loci') + xlab('QTL') + 
  geom_text(data=AD_total_by_qtl, aes(x=QTL, y=total, label=total, fill = NULL, alpha = NULL, size=15), nudge_y = 1) + 
  facet_grid(. ~celltype, scales = "free_x", space = "free_x") + 
  theme(strip.text.x = element_text(size = 17), 
        axis.title = element_text(size=17, colour = "black"), a
        xis.text.y = element_text(size=15, colour="black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size=15), 
        legend.spacing = unit(1, 'cm'), 
        legend.text = element_text(size=15, colour = "black"), 
        legend.title = element_text(size=15, colour = "black")) + 
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme(strip.placement = "outside", 
        strip.background = element_blank(),
        strip.text =  element_text(size = 10)) +
  ggtitle(gwas_name)
show(barplot)

ggsave(filename, plot = barplot, 
       width=30 ,height=12, 
       dpi=600, limitsize = FALSE)


## COLOC dot plot

AD_total_colocs$PP.H4.abf <- as.numeric(AD_total_colocs$PP.H4.abf)
AD_total_colocs$PP.H4.abf <- AD_total_colocs$PP.H4.abf <- round(AD_total_colocs$PP.H4.abf, digits=2)
coloc_combined_melted <- melt(AD_total_colocs, ID='locus')
coloc_combined_melted_filt <-  coloc_combined_melted[order(coloc_combined_melted[,'new_QTL_Gene'],coloc_combined_melted[,'type'], -coloc_combined_melted[,'value']),]
df <- distinct(coloc_combined_melted_filt, locus, new_QTL_Gene, type, .keep_all= TRUE)
df$celltype <- factor(df$celltype , levels = dataset_celltype)

subsetted <- df

p1 <- ggplot(data = subsetted, aes(x=QTL, y=new_QTL_Gene, color=celltype)) +
  geom_point(aes(size = value, alpha = value), position = position_nudge(x = -0.2)) +
  scale_colour_manual(values=colorset) +
  theme_minimal() +
  geom_text(aes(label=value),color = 'black',nudge_x = 0.2, size = 4) +
  facet_grid(locus ~ ., scales = "free_y",space = "free_y" )+
  ylab('') +
  xlab('') +
  theme(strip.text.y = element_text(angle = 0, face='italic', colour = "black", hjust = 0, size=13),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size=13),
        axis.text.y = element_text(face = "bold.italic", colour = "black",size=13),
        strip.background = element_blank(),
        legend.position = "right",
        legend.box.margin = margin(c(0,0,0,0)),
        legend.margin = margin(c(0,0,0,0)),
        strip.placement = "outside",
        panel.spacing.y = unit(x = 0,units = "points"), 
        panel.border = element_rect(fill = NA, size = 0.1),
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black")) +
  facet_grid(locus ~ ., scales = "free_y",space = "free_y" ) +
  ggtitle(gwas_name)

p1 <- p1 + facet_grid(locus ~ celltype, scales = "free", space = "free") + 
  theme(strip.text = element_text(size = 13, color = "black"))

show(p1)

ggsave(filename, plot = p1, 
       width=25 ,height=40, 
       dpi=600)



