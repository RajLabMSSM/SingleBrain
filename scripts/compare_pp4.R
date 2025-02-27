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


# COLOC compare 
COLOC_res <- read_tsv('~/all_COLOC_results_merged_H4_0.5_with_LD_filtered.tsv.gz')

AD_total_colocs <- COLOC_res

genes <- AD_total_colocs$QTL_Gene[grepl('ENSG0', AD_total_colocs$QTL_Gene, fixed=TRUE)]
genes <- sapply(strsplit(genes, '\\.'), '[[', 1)
AD_total_colocs$QTL_Ensembl[is.na(AD_total_colocs$QTL_Ensembl)] = AD_total_colocs$QTL_Gene[is.na(AD_total_colocs$QTL_Ensembl)]
AD_total_colocs$QTL_Ensembl <- sapply(strsplit(AD_total_colocs$QTL_Ensembl, '\\.'), '[[', 1)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
colnames(geneIDs1) = c('new_QTL_Gene', 'QTL_Ensembl')
AD_total_colocs <- left_join(AD_total_colocs, geneIDs1, by='QTL_Ensembl')
AD_total_colocs$new_QTL_Gene[is.na(AD_total_colocs$new_QTL_Gene)] = AD_total_colocs$QTL_Gene[is.na(AD_total_colocs$new_QTL_Gene)]
AD_total_colocs <- AD_total_colocs[, c('disease',"GWAS","locus", 'GWAS_P',"PP.H4.abf","new_QTL_Gene","QTL", 'celltype')]

AD_total_colocs <- AD_total_colocs[which(AD_total_colocs$QTL %in% snrna_qtl ),]
AD_total_colocs <- AD_total_colocs[!is.na(AD_total_colocs$PP.H4.abf),]

AD_total_colocs <- AD_total_colocs %>% group_by(GWAS, locus, new_QTL_Gene, QTL) %>% slice_max(PP.H4.abf)


COLOC_qtl1 <- AD_total_colocs %>% dplyr::filter((QTL==qtl1_name))
COLOC_qtl2 <- AD_total_colocs %>% dplyr::filter((QTL==qtl2_name))

COLOC_qtl1 <- COLOC_qtl1[which(COLOC_qtl1$new_QTL_Gene %in% COLOC_qtl2$new_QTL_Gene),]
COLOC_qtl1 <- COLOC_qtl1[which(COLOC_qtl1$locus %in% COLOC_qtl2$locus),]
COLOC_qtl1 <- COLOC_qtl1 %>% dplyr::select(-GWAS, -GWAS_P)
COLOC_qtl2 <- COLOC_qtl2[which(COLOC_qtl2$new_QTL_Gene %in% COLOC_qtl1$new_QTL_Gene),]
COLOC_qtl2 <- COLOC_qtl2[which(COLOC_qtl2$locus %in% COLOC_qtl1$locus),]
COLOC_qtl2 <- COLOC_qtl2 %>% dplyr::select(-GWAS, -GWAS_P)

COLOC_qtl <- inner_join(COLOC_qtl1, COLOC_qtl2, by=c('GWAS', "disease",'locus',
                                                     'new_QTL_Gene', 'celltype'))
COLOC_qtl <- COLOC_qtl%>% group_by(GWAS, new_QTL_Gene) %>% slice_max(PP.H4.abf.x)
COLOC_qtl <- COLOC_qtl%>% group_by(GWAS, new_QTL_Gene) %>% slice_max(PP.H4.abf.y)

options(ggrepel.max.overlaps = Inf)
COLOC_qtl_fill <- COLOC_qtl%>% dplyr::filter(((PP.H4.abf.x>0.5)|(PP.H4.abf.y>0.5)))
COLOC_qtl_fill$ratio <-  log2(COLOC_qtl_fill$PP.H4.abf.x /COLOC_qtl_fill$PP.H4.abf.y)
COLOC_qtl_fill <- COLOC_qtl_fill %>% mutate(ratio_color = case_when(ratio <= -5 ~ -5,
                                                                    ratio >= 5 ~ 5,
                                                                    (ratio > -5)&(ratio < 5) ~ ratio))

COLOC_qtl_fill_tmp <- COLOC_qtl_fill %>% dplyr::filter(disease==dis)
COLOC_qtl_diff <- COLOC_qtl_fill_tmp%>% dplyr::filter(((ratio_color < -1)|
                                                         ((PP.H4.abf.x <= 0.8)&(PP.H4.abf.y >= 0.5))|
                                                         ((PP.H4.abf.x >= 0.8)|(PP.H4.abf.y <= 0.5))
))

p1 <- ggplot(COLOC_qtl_fill_tmp, aes(x=PP.H4.abf.x, y=PP.H4.abf.y, 
                                     color=ratio_color, label=new_QTL_Gene)) +
  geom_point(size=4) + 
  scale_color_continuous(type = "viridis", 
  ) +
  labs(title = dis,
       x=paste0(qtl1_name,' PP.H4') , y = 'MiGA PP.H4'
  ) +
  scale_x_continuous(limits = c(0.0,1.05),
                     breaks = c(0,0.5,1)) +
  scale_y_continuous(limits = c(0.0,1.05),
                     breaks = c(0,0.5,1)) +
  geom_vline(xintercept = c(0.9,1), linetype="dashed", 
             color = "black", size=0.5, alpha=0.5 ) +
  geom_hline(yintercept = c(0.9,1), linetype="dashed", 
             color = "black", size=0.5, alpha=0.5 ) +
  theme_classic() +
  theme( axis.text = element_text(size=12, color='black'),
         axis.title = element_text(size=12, color='black', face='bold'),
         legend.position = 'right',
         legend.text = element_blank(),
         legend.title = element_blank(),
         legend.ticks = element_blank(),
         strip.background = element_blank(), 
         strip.text = element_text(size=12, color='black', face='bold'),
  ) 
p1 <- p1 + ggrepel::geom_text_repel(data = COLOC_qtl_fill_tmp, aes(label=new_QTL_Gene),
                                    color='black', size=3)

show(p1)

ggsave(filename, plot = p1, 
       width=5.5 ,height=5, 
       dpi=600)



