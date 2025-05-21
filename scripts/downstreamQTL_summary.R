
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



# COLOC summary
COLOC_total <- read_tsv(paste0(path_dir,"/all_COLOC_results_merged_H4_0.5_with_LD.tsv.gz"))

# COLOC summary stack bar
AD_total_colocs <- COLOC_result

AD_total_colocs <- AD_total_colocs %>% 
  dplyr::filter(grepl('SingleBrain_', AD_total_colocs$QTL)) %>% 
  dplyr::filter(PP.H4.abf > 0.8) 

snrna_qtl <- c('SingleBrain_Ast','SingleBrain_MG','SingleBrain_Ext', 'SingleBrain_IN',
               'SingleBrain_OD', 'SingleBrain_OPC','SingleBrain_End'
)
gwas_list <- c('VanRheenenEUR_2021','Daner_2020','TrubetskoyEUR_2022',
               'Nalls23andMe_2019','Bellenguez_2021', 'IMSGC_2019')

### Extract the GWAS and QTL to plot
AD_total_colocs <- AD_total_colocs[which(AD_total_colocs$QTL %in% snrna_qtl ),]
AD_total_colocs <- AD_total_colocs[which(AD_total_colocs$GWAS %in% gwas_list ),]


AD_total_colocs_summary <- AD_total_colocs %>% 
  arrange(QTL, desc(PP.H4.abf)) %>% 
  group_by(QTL, celltype) %>% 
  distinct(locus, .keep_all = TRUE) %>% count(GWAS)


AD_total_colocs_summary$GWAS <- gsub('Bellenguez_2021','Bellenguez_2022',AD_total_colocs_summary$GWAS)
AD_total_colocs_summary$GWAS <- gsub('TrubetskoyEUR_2022','Trubetskoy2022',AD_total_colocs_summary$GWAS)
AD_total_colocs_summary$GWAS <- gsub('Daner_2020','Mullins_2021',AD_total_colocs_summary$GWAS)
AD_total_colocs_summary$GWAS <- gsub('VanRheenenEUR_2021','Rheenen_2021',AD_total_colocs_summary$GWAS)


gwas_list <- c('Yengo_2022','Ishigaki_2022','Rheenen_2021',
               'IMSGC_2019','Mullins_2021','Trubetskoy2022',
               'Nalls23andMe_2019','Marioni_2018','Jansen_2018', 
               'Kunkle_2019','Bellenguez_2022')

disease_list <- c('ALS','MS', 'BP','SCZ','PD','AD')

cell_list <- c('Endothelial cells','OPCs','Microglia',
               'Astrocytes','Oligodendrocytes','Inhibitory Neurons', 'Excitatory Neurons'
)

### Annotate the disease with GWAS name
AD_total_colocs_summary[AD_total_colocs_summary$GWAS=='Bellenguez_2022','Disease'] <- 'AD'
AD_total_colocs_summary[AD_total_colocs_summary$GWAS=='Kunkle_2019','Disease'] <- 'AD'
AD_total_colocs_summary[AD_total_colocs_summary$GWAS=='Jansen_2018','Disease'] <- 'AD'
AD_total_colocs_summary[AD_total_colocs_summary$GWAS=='Marioni_2018','Disease'] <- 'AD'
AD_total_colocs_summary[AD_total_colocs_summary$GWAS=='Nalls23andMe_2019','Disease'] <- 'PD'
AD_total_colocs_summary[AD_total_colocs_summary$GWAS=='Trubetskoy2022','Disease'] <- 'SCZ'
AD_total_colocs_summary[AD_total_colocs_summary$GWAS=='Mullins_2021','Disease'] <- 'BP'
AD_total_colocs_summary[AD_total_colocs_summary$GWAS=='Rheenen_2021','Disease'] <- 'ALS'
AD_total_colocs_summary[AD_total_colocs_summary$GWAS=='IMSGC_2019','Disease'] <- 'MS'
AD_total_colocs_summary[AD_total_colocs_summary$GWAS=='Ishigaki_2022','Disease'] <- 'RA'
AD_total_colocs_summary[AD_total_colocs_summary$GWAS=='Yengo_2022','Disease'] <- 'Height'

AD_total_colocs_summary$GWAS <- factor(AD_total_colocs_summary$GWAS , levels = gwas_list)
AD_total_colocs_summary$celltype <- factor(AD_total_colocs_summary$celltype , levels = cell_list)
AD_total_colocs_summary$Disease <- factor(AD_total_colocs_summary$Disease , levels = disease_list)

AD_total_colocs_summary <- AD_total_colocs_summary %>%
  group_by(GWAS) %>%
  mutate(Freq = n/sum(n))

AD_total_colocs_summary$Freq <- AD_total_colocs_summary$Freq * 100
AD_total_colocs_summary$Freq  <- round(AD_total_colocs_summary$Freq, digits=2)

### plotting
p_stack <-ggplot(AD_total_colocs_summary, aes(fill=celltype, y=n, x=Disease, label = celltype)) + 
  geom_bar(position="fill", stat="identity", color = "black") + 
  scale_fill_manual(values = colorset) + 
  scale_x_discrete(position = "top") + 
  theme_minimal() +
  ylab('Celltype Proportion') +
  theme(axis.text.y = element_text(angle = 0, color = "black", 
                                   face='italic',hjust = 1, 
                                   colour = "black", size=13),
    strip.placement = "outside",
    axis.title.x = element_text(vjust=0)
  ) + 
  coord_flip()

show(p_stack)

ggsave(filename, plot = p_stack, 
       width=8 ,height=5, 
       dpi=600)


## COLOC MR summary
COLOC_total <- read_tsv(paste0(path_dir,"/all_COLOC_results_merged_H4_0.5_with_LD.tsv.gz"))

### Load MR result and filter by significant value
MR_total <- read_tsv(paste0(path_dir,"/all_MR_results_merged_H4_1_snp_wFDR.tsv.gz"))

MR_total <- MR_total %>% dplyr::filter(grepl("Inverse variance weighted", method)) 

MR_result <- MR_total %>% group_by(GWAS) %>% 
  mutate(FDR = p.adjust (pval, method='BH'),
         bonf = p.adjust(pval, method = "bonferroni")
  )

MR_fdr <- MR_result %>% dplyr::filter(FDR < 0.05)

snrna_qtl <- c('SingleBrain_MG',
               'SingleBrain_Ast',
               'SingleBrain_Ext',
               'SingleBrain_IN',
               'SingleBrain_OD',
               'SingleBrain_OPC',
               'SingleBrain_End'
               # 'MiGA'
)
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

AD_total_mrs$QTL <- gsub('SingleBrain_','',AD_total_mrs$QTL)

### Extract necessary columns
subsetted2 <- AD_total_mrs %>% dplyr::select(GWAS, locus, new_QTL_Gene, QTL, celltype, value)
subsetted2$method <- 'MR_positive'
subsetted2[(subsetted2$value<0),'method'] = 'MR_negative'
subsetted2$value <- 1
subsetted2$type <- paste0('MR_', subsetted2$QTL)

subsetted <- COLOC_result %>% dplyr::select(GWAS, locus, new_QTL_Gene, QTL, celltype, value=PP.H4.abf)
subsetted$type <- paste0('COLOC_', subsetted$QTL)
subsetted$method <- 'COLOC'
subsetted4 <- rbind(subsetted2[which(subsetted2$id %in% subsetted$id ),],subsetted)

### Make MR sig or not
subsetted5 <- subsetted %>% dplyr::filter(value>=0.8)
subsetted5$MR <- 'No'
subsetted5[which(subsetted5$id %in% subsetted2$id ),'MR'] = 'Yes'
table(subsetted5$MR, subsetted5$GWAS)

GWAS_loci <- data.frame(Bellenguez_2021 = c(75),
                        Daner_2020 = c(64),
                        IMSGC_2019 = c(137),
                        Nalls23andMe_2019 = c(71),
                        TrubetskoyEUR_2022 = c(217),
                        VanRheenenEUR_2021 = c(15)
)

GWAS_loci <- GWAS_loci %>% t() %>% as.data.frame()
colnames(GWAS_loci) <- c('counts')
GWAS_loci$GWAS <- rownames(GWAS_loci)
rownames(GWAS_loci) <- NULL
GWAS_loci$type <- 'Total'

subsetted6 <- subsetted5 %>% arrange(desc(MR)) %>% 
  dplyr::select(GWAS, locus, MR) %>% 
  distinct(GWAS, locus,.keep_all=TRUE)
subsetted6

df_count <- table(subsetted6$MR, subsetted6$GWAS) %>% 
  t() %>% as.matrix() %>% 
  as.data.frame()
colnames(df_count) <- c('GWAS','type','counts')

df_count_coloc <- subsetted6 %>% 
  distinct(GWAS, locus) %>% 
  group_by(GWAS)  %>% 
  summarise(counts=n())
df_count_coloc$type <- 'coloc'

df_count <- rbind(GWAS_loci, df_count)
df_count <- rbind(df_count_coloc, df_count)

df_count$disease <- 'AD'
df_count[df_count$GWAS=='Daner_2020', 'disease'] =  'BPD'
df_count[df_count$GWAS=='IMSGC_2019', 'disease'] =  'MS'
df_count[df_count$GWAS=='Nalls23andMe_2019', 'disease'] =  'PD'
df_count[df_count$GWAS=='TrubetskoyEUR_2022', 'disease'] =  'SCZ'
df_count[df_count$GWAS=='VanRheenenEUR_2021', 'disease'] =  'ALS'

### Categorize the result
df_count_seq <- df_count %>% dplyr::filter(type=='Total') %>% arrange(desc(counts))

df_count$disease <- factor(df_count$disease, levels = df_count_seq$disease)

df_count_total <- df_count %>% dplyr::filter(type=='Total') %>% dplyr::select(GWAS, counts)
colnames(df_count_total) <- c('GWAS','Total')
df_count_coloc <- df_count %>% dplyr::filter(type=='coloc') %>% dplyr::select(GWAS, counts)
colnames(df_count_coloc) <- c('GWAS','COLOC')
df_count_mr <- df_count %>% dplyr::filter(type=='Yes') %>% dplyr::select(GWAS, counts)
colnames(df_count_mr) <- c('GWAS','MR')

df_count_label <- merge(df_count_total, df_count_coloc)
df_count_label <- merge(df_count_label, df_count_mr)
df_count_label$label <- paste(df_count_label$MR, df_count_label$COLOC, df_count_label$Total) 
df_count_label <- df_count_label %>% dplyr::select(GWAS, label)
df_count <- merge(df_count, df_count_label, all.x=TRUE)

### plotting
p_overlap_mr <- df_count %>% dplyr::filter(type=='Total') %>% 
  ggplot(aes(x=disease, y=counts)) + 
  geom_bar(stat = "identity", fill='darkgrey') +
  geom_bar(data=df_count %>% dplyr::filter(type=='coloc'), 
           aes(x=disease, y=counts), 
           fill='#cab2d6',
           stat = "identity") +
  geom_text(aes(label=label),
            position = position_dodge(width=-1),  
            size=5, hjust = -0.1, color='black'
  ) +
  geom_bar(data=df_count %>% dplyr::filter(type=='Yes'), 
           aes(x=disease, y=counts), 
           fill='#6a3d9a',
           stat = "identity") +
  scale_y_continuous(limits = c(0, 260), ) +
  theme_bw() +
  coord_flip() +
  labs(title='GWAS discovery with MR',
       x='Disease', y='GWAS loci') +
  theme(axis.title = element_text(size = 12, color='black'),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.text = element_text(size = 12, color='black'),
        plot.title = element_text(size = 12, color='black', face='bold')
  )

p_overlap_mr

ggsave(filename, plot = p_overlap_mr, 
       width=7 ,height=5, 
       dpi=600)


## COLOC fine-mapping
### Load the fine-mapping result
finemap_all_res <- read_tsv('~/all_results_merged_PIP_0.1.tsv.gz') 
finemap_all_res_sim_095 <- finemap_all_res[finemap_all_res$PIP >= 0.95, 
                                           c('feature', 'QTL')] %>% unique()
finemap_all_res_sim_095$PIP095 <- 'Yes'

### Extract coloc result
subsetted7 <- subsetted %>% dplyr::filter(value>=0.8)
subsetted7 <- left_join(subsetted7, tmp_res_095, by=c('new_QTL_Gene'='Symbol','QTL'))
subsetted7[is.na(subsetted7)]='No'

subsetted8 <- subsetted7 %>% arrange(desc(PIP095)) %>% 
  dplyr::select(GWAS, locus, PIP095) %>% 
  distinct(GWAS, locus,.keep_all=TRUE)

df_count2 <- table(subsetted8$PIP095, subsetted8$GWAS) %>% 
  t() %>% as.matrix() %>% as.data.frame()
colnames(df_count2) <- c('GWAS','type','counts')

df_count_coloc <- subsetted8 %>% distinct(GWAS, locus) %>% 
  group_by(GWAS)  %>% summarise(counts=n())
df_count_coloc$type <- 'coloc'

df_count2 <- rbind(GWAS_loci, df_count2)
df_count2 <- rbind(df_count_coloc, df_count2)

df_count2$disease <- 'AD'
df_count2[df_count2$GWAS=='Daner_2020', 'disease'] =  'BPD'
df_count2[df_count2$GWAS=='IMSGC_2019', 'disease'] =  'MS'
df_count2[df_count2$GWAS=='Nalls23andMe_2019', 'disease'] =  'PD'
df_count2[df_count2$GWAS=='TrubetskoyEUR_2022', 'disease'] =  'SCZ'
df_count2[df_count2$GWAS=='VanRheenenEUR_2021', 'disease'] =  'ALS'

### Categorize the result
df_count_seq <- df_count2 %>% dplyr::filter(type=='Total') %>% arrange(desc(counts))

df_count2$disease <- factor(df_count2$disease, levels = df_count_seq$disease)

df_count_total <- df_count2 %>% dplyr::filter(type=='Total') %>% dplyr::select(GWAS, counts)
colnames(df_count_total) <- c('GWAS','Total')
df_count_coloc <- df_count2 %>% dplyr::filter(type=='coloc') %>% dplyr::select(GWAS, counts)
colnames(df_count_coloc) <- c('GWAS','COLOC')
df_count_mr <- df_count2 %>% dplyr::filter(type=='Yes') %>% dplyr::select(GWAS, counts)
colnames(df_count_mr) <- c('GWAS','Finemap')

df_count_label <- merge(df_count_total, df_count_coloc)
df_count_label <- merge(df_count_label, df_count_mr)

df_count_label$label <- paste(df_count_label$Finemap, df_count_label$COLOC, df_count_label$Total) 

df_count_label <- df_count_label %>% dplyr::select(GWAS, label)
df_count2 <- merge(df_count2, df_count_label, all.x=TRUE)

### plotting
p_overlap_finemap <- df_count2 %>% dplyr::filter(type=='Total') %>% ggplot(aes(x=disease, y=counts)) + 
  geom_bar(stat = "identity", fill='darkgrey') +
  geom_bar(data=df_count2 %>% dplyr::filter(type=='coloc'), 
           aes(x=disease, y=counts), 
           fill='#a6cee3',
           stat = "identity") +
  geom_text(aes(label=label),
            position = position_dodge(width=-1),  
            size=5, hjust = -0.1, color='black'
  ) +
  geom_bar(data=df_count2 %>% dplyr::filter(type=='Yes'), 
           aes(x=disease, y=counts), 
           fill='#1f78b4',
           stat = "identity") +
  scale_y_continuous(limits = c(0, 260), ) +
  theme_classic() +
  coord_flip() +
  labs(title='GWAS discovery with Fine-mapping (PIP > 0.95)',
       x='', y='Number of GWAS loci') +
  theme(axis.title = element_text(size = 12, color='black'),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.text = element_text(size = 12, color='black'),
        plot.title = element_text(size = 12, color='black', face='bold')
  )

p_overlap_finemap

ggsave(filename, plot = p_overlap_finemap, 
       width=8 ,height=5,
       dpi=600)














