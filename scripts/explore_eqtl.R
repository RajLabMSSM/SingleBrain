library(readr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(ggplotify)
library(data.table)

# Loading eGene
dirs <- list.dirs(here::here(paste0(path_mmqtl,"mmQTL/")),full.names = FALSE, recursive = FALSE)
assoc <- here::here(paste0(path_mmqtl, "mmQTL/",dirs,"/", dirs, "_top_assoc.tsv.gz"))
meta <- here::here(paste0(path_mmqtl, "mmQTL/",dirs,"/","phenotype_metadata.tsv"))

file.exists(assoc)
file.exists(meta)

res <- map2( assoc, meta, ~{
  a <- read_tsv(.x)
  m <- read_tsv(.y, col_names = c("chr", "start", "end", "feature"))
  d <- inner_join(m, a, by = c("feature", "chr"))
})

names(res) <- dirs

qtl_summary <- map_df( res, ~{
  sig <- filter(.x , qval < 0.05)
  tibble( 
    features = nrow(.x), 
    sig_features = nrow(sig),
    null_features = nrow(.x) - nrow(sig),
    perc_features = signif(nrow(sig) / nrow(.x) *  100, 2),
    # genes = length(unique(.x$group)),
    # sig_genes = length(unique(sig$group)),
    # perc_genes = signif(length(unique(sig$group)) / length(unique(.x$group)) * 100, 2)
  )
}, .id = "dataset") %>%
  tidyr::separate(dataset, into = c("reference", "phenotype"), sep = "_", extra = "merge") 


# eGene by cell proportion
egene_number_pro_plot <- ggplot(data = egenes_number_sb, 
                                aes(x = Freq, y = `eGENEs`, 
                                    color = as.factor(cell_type))) + 
  geom_point(stat = "identity", size = 5.5) +
  scale_color_manual(values = colorset) +
  theme_bw() + 
  scale_y_continuous(name = "eGenes (n)", breaks = c(4000,8000,10000,12000), limits = c(0,12500)) +
  scale_x_continuous(name = "Cell proportion (%)", breaks = c(1, 5, 10,15,25, 40, 50)) +
  theme(axis.text.x = element_text( hjust = 0.5, vjust = 0.9, size = 13, color='black'),
        axis.text.y = element_text(angle = 45, size = 13, color='black'),
        legend.position = 'none'
  ) +
  labs(title = "Number of significant eGENEs found", 
       x = "Cell proportion", 
       y = "# eGENEs", 
       color = "Cell Type\n") +
  labs(color = "") +
  ggrepel::geom_text_repel(data = egenes_number_sb, aes(label=celltype),
                           force=1, point.padding = unit(0.5,'lines'),
                           box.padding=unit(0.5,'lines'),
                           direction = "y", hjust = 1, vjust = 1, 
                           color='black' )  + 
  geom_smooth(method = lm, formula = eGENEs~Freq, data = egenes_number_sb)
show(egene_number_pro_plot) 

ggsave(filename, plot = egene_number_pro_plot, 
       width=5 ,height=5, dpi=600, 
       limitsize = FALSE)


# eGene bar
egenes_number_plot <- egenes_number %>% dplyr::filter(cohort=='SingleBrain')
egenes_number_plot$label <- gsub("SingleBrain_","",egenes_number_plot$label )
p <- ggplot(data=egenes_number_plot %>% dplyr::filter(cohort=='SingleBrain'), 
            aes(x=reorder(label, -eGENEs), y=eGENEs, fill=factor(cell_type))) +
  geom_bar(stat="identity") +
  theme_classic() +
  xlab("Cell type")  + 
  ylab("eGenes")  + 
  theme(legend.position = "none") +
  ggtitle("SingleBrain eGene by Celltype ") + theme(plot.title = element_text(size = 14, hjust = 0.5, face='bold')) +
  theme(axis.text.x = element_text( size = 12, color='black', angle = 45, vjust=1, hjust=1),
        axis.text.y = element_text( size = 12,  color='black', vjust=1, hjust=1) ,
        axis.title = element_text( size = 12,  color='black'),
  ) +
  geom_text( aes(label = eGENEs), vjust = -1, size = 4, face='bold') +
  scale_fill_manual(values=colorset)
p

ggsave(filename, plot = p, width=8 ,height=7, 
       dpi=600, limitsize = FALSE)


# SingleBrain comparison
egene_number_plot <- ggplot(data = egenes_number, 
                            aes(x = cohort, y = `eGENEs`, 
                                color = as.factor(cell_type), 
                                shape=cohort)) + 
  geom_point(stat = "identity", size = 5.5) +
  scale_shape_manual(values = c(17, 15, 16))+
  scale_color_manual(values = colorset) +
  theme_bw() + 
  scale_y_continuous(name = "eGenes (N)", breaks = c(100, 2000,4000,8000,12000)) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.9, size = 13, color='black'),
        axis.text.y = element_text(angle = 45, size = 13, color='black'),
        axis.title = element_text(size = 13, color='black'),
  ) + 
  labs(title = "Number of significant eGENEs found", 
       x = "Cohorts", y = "# eGENEs", color = "Cell Type\n") +
  labs(color = "")
show(egene_number_plot)

ggsave(filename, plot = egene_number_plot, 
       width=6.5 ,height=5, 
       dpi=600, limitsize = FALSE)


# upset plot

mmqtl_ast_sig <- mmqtl_ast_top %>% mutate(Ast=case_when(qval < 0.05 ~ 1, qval >= 0.05 ~ 0)) %>% dplyr::select(feature, Ast)
mmqtl_end_sig <- mmqtl_end_top %>% mutate(End=case_when(qval < 0.05 ~ 1, qval >= 0.05 ~ 0)) %>% dplyr::select(feature, End) 
mmqtl_ext_sig <- mmqtl_ext_top %>% mutate(Ext=case_when(qval < 0.05 ~ 1, qval >= 0.05 ~ 0)) %>% dplyr::select(feature, Ext)
mmqtl_in_sig <- mmqtl_in_top %>% mutate(IN=case_when(qval < 0.05 ~ 1, qval >= 0.05 ~ 0)) %>% dplyr::select(feature, IN)
mmqtl_mg_sig <- mmqtl_mg_top %>% mutate(MG=case_when(qval < 0.05 ~ 1, qval >= 0.05 ~ 0)) %>% dplyr::select(feature, MG)
mmqtl_od_sig <- mmqtl_od_top %>% mutate(OD=case_when(qval < 0.05 ~ 1, qval >= 0.05 ~ 0)) %>% dplyr::select(feature, OD)
mmqtl_opc_sig <- mmqtl_opc_top %>% mutate(OPC=case_when(qval < 0.05 ~ 1, qval >= 0.05 ~ 0)) %>% dplyr::select(feature, OPC)


mmqtl_upset <- merge(mmqtl_ast_sig, mmqtl_end_sig, all=TRUE, by = 'feature')
mmqtl_upset <- merge(mmqtl_upset, mmqtl_ext_sig, all=TRUE, by = 'feature')
mmqtl_upset <- merge(mmqtl_upset, mmqtl_in_sig, all=TRUE, by = 'feature')
mmqtl_upset <- merge(mmqtl_upset, mmqtl_mg_sig, all=TRUE, by = 'feature')
mmqtl_upset <- merge(mmqtl_upset, mmqtl_od_sig, all=TRUE, by = 'feature')
mmqtl_upset <- merge(mmqtl_upset, mmqtl_opc_sig, all=TRUE, by = 'feature')
rownames(mmqtl_upset) <- mmqtl_upset$feature
mmqtl_upset <- mmqtl_upset  %>% dplyr::select(-feature) 

mmqtl_upset[is.na(mmqtl_upset)] <- 0

mmqtl_upset_total_count <- colSums(mmqtl_upset) %>% as.data.frame()
mmqtl_upset_total_count$celltype <- rownames(mmqtl_upset_total_count)
rownames(mmqtl_upset_total_count) <- NULL
colnames(mmqtl_upset_total_count) <- c('overlap','celltype')
mmqtl_upset_total_count <- mmqtl_upset_total_count %>% arrange(desc(overlap))
mmqtl_upset_total_count$colors <- mmqtl_upset_total_count$celltype

mmqtl_upset_celltype <- mmqtl_upset[rowSums(mmqtl_upset) ==1,]

mmqtl_upset_celltype_count <- colSums(mmqtl_upset_celltype) %>% as.data.frame()
mmqtl_upset_celltype_count$celltype <- rownames(mmqtl_upset_celltype_count)
rownames(mmqtl_upset_celltype_count) <- NULL
colnames(mmqtl_upset_celltype_count) <- c('uniq','celltype')

mmqtl_upset_total_count <- left_join(mmqtl_upset_total_count, mmqtl_upset_celltype_count, by='celltype')

cell_list <- c('Ast','End','Ext','IN','MG','OD','OPC')

df_sec <- data.frame()

for(num1 in 1:length(cell_list)){
  for(num2 in 1:length(cell_list)){
    if(num1>=num2){next}
    celltype1 <- cell_list[num1]
    celltype2 <- cell_list[num2]
    term <- paste0(celltype1,'-',celltype2)
    df_upset_tmp <- mmqtl_upset %>% dplyr::select(celltype1, celltype2)
    df_upset_tmp_celltype <- df_upset_tmp[rowSums(df_upset_tmp) ==2,]
    df_sec_tmp <- data.frame(celltype=c(term), 
                             overlap=c(nrow(df_upset_tmp_celltype))
    )
    df_sec <- rbind(df_sec, df_sec_tmp)
  }
}
df_sec <- df_sec %>% arrange(desc(overlap))
df_sec$colors <- 'Second'
df_sec$uniq <- 0

df_list <- c(mmqtl_upset_total_count$celltype, df_sec$celltype)

df_result <- rbind(mmqtl_upset_total_count, df_sec)
df_result$uniqbar <- df_result$overlap - df_result$uniq

df_result$celltype <- factor(df_result$celltype, level=df_list)

### Plotting
p_bar <- df_result %>% ggplot(aes(x=celltype, y=overlap)) + 
  geom_bar( width = 0.8, fill='darkgrey', stat="identity") +
  geom_text(aes(label=overlap), 
            hjust=-0.15, angle=90, size = 4, color='black') +
  geom_bar(data=df_result, aes(x=celltype, y=uniqbar, fill=colors), width = 0.8, stat="identity") + 
  scale_y_continuous(limits = c(0,13000)) +
  scale_fill_manual(values = color_upset) +
  theme_minimal() +
  xlab('') + ylab('') +
  theme(axis.title = element_text(size=12, color='black'),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = 'none'
  )

p_bar

### Making total bar
df_result$point <- df_result$celltype
df_result_point <- separate_rows(df_result, point, sep = '-')

df_point_back <- data.frame()
for(celltype in df_result$celltype){
  df_point_back_tmp <- data.frame(celltype=rep(celltype, 7),
                                  term=cell_list
  )
  df_point_back <- rbind(df_point_back, df_point_back_tmp)
}

df_point_back$term <- factor(df_point_back$term, 
                             levels = mmqtl_upset_total_count$celltype)
df_result_point$point <- factor(df_result_point$point, 
                                levels = mmqtl_upset_total_count$celltype)

p_point <-
  df_point_back %>% ggplot(aes(x = celltype, y = term))  +
  geom_point(size = 6.5, color='#d4d4d4') +
  geom_point(data=df_result_point, 
             aes(x = celltype, y = point, color=point), 
             size = 6.5) +
  scale_color_manual(values = color_upset) +
  theme_void() +
  xlab('') + ylab('') +
  theme(axis.title = element_text(size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.text.x = element_blank(),
        legend.position = 'none'
  ) 

p_point  

p_com <- cowplot::plot_grid(p_bar, p_point,
                            align='v', nrow=2, rel_heights = c(4,  2))
p_com

ggsave(filename, plot = p_com, width=8 ,height=6, dpi=600, limitsize = FALSE)


# Subtypes

## eGene with gray celltype specific
dirs <- list.dirs(here::here(paste0(path_mmqtl,"mmQTL/")),full.names = FALSE, recursive = FALSE)
assoc <- here::here(paste0(path_mmqtl, "mmQTL/",dirs,"/", dirs, "_top_assoc.tsv.gz"))
meta <- here::here(paste0(path_mmqtl, "mmQTL/",dirs,"/","phenotype_metadata.tsv"))

file.exists(assoc)
file.exists(meta)

res <- map2( assoc, meta, ~{
  a <- read_tsv(.x)
  m <- read_tsv(.y, col_names = c("chr", "start", "end", "feature"))
  d <- inner_join(m, a, by = c("feature", "chr"))
})

names(res) <- dirs

qtl_summary <- map_df( res, ~{
  sig <- filter(.x , qval < 0.05)
  tibble( 
    features = nrow(.x), 
    sig_features = nrow(sig),
    null_features = nrow(.x) - nrow(sig),
    perc_features = signif(nrow(sig) / nrow(.x) *  100, 2),
    # genes = length(unique(.x$group)),
    # sig_genes = length(unique(sig$group)),
    # perc_genes = signif(length(unique(sig$group)) / length(unique(.x$group)) * 100, 2)
  )
}, .id = "dataset") %>%
  tidyr::separate(dataset, into = c("reference", "phenotype"), sep = "_", extra = "merge") 


p <- ggplot(data=egenes_number_sb_only , aes(x=reorder(celltype, -eGENEs), y=eGENEs, 
)) +
  geom_bar(stat="identity", fill='grey') +
  geom_bar(data=egenes_number_sb_only , aes(x=reorder(celltype, -major_eGENEs), y=major_eGENEs, 
                                            fill=factor(cell_type)), stat="identity") +
  geom_text(aes(label=eGENEs), 
            hjust=0, angle=90, size = 4, color='black') +
  facet_grid(. ~cell_type, scales = "free_x", space = "free_x") +
  theme_classic() +
  xlab("Cell type")  + 
  ylab("eGenes")  + 
  scale_y_continuous(limits = c(0,13000)) +
  ggtitle("") + 
  theme(legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0.5, face='bold'),
        axis.text.x = element_text( size = 12, angle = 45, vjust=1, hjust=1, color='black'),
        axis.text.y = element_text( size = 12,  vjust=1, hjust=1, color='black'),
        axis.title = element_text( size = 12,  color='black'),
        strip.background = element_blank(),
        strip.text = element_blank()
  ) +
  scale_fill_manual(values=colorset)
p

ggsave(filename, plot = p, width=10 ,height=4.5, 
       dpi=600, limitsize = FALSE)

## celltype specific only
egenes_number_sb_only$cell_type <- factor(egenes_number_sb_only$cell_type , levels = cellorder)

p <- ggplot(data=egenes_number_sb_only , aes(x=reorder(celltype, -eGENEs), y=sub_only_counts, 
                                             fill=factor(cell_type))) +
  geom_bar(stat="identity") +
  geom_text(aes(label=sub_only_counts),
            hjust=0, angle=90, size = 4, color='black') +
  facet_grid(. ~cell_type, scales = "free_x", space = "free_x") +
  theme_classic() +
  xlab("Cell type")  + 
  ylab("eGenes")  + 
  scale_y_continuous(limits = c(0,3200)) +
  ggtitle("mmQTL eGene by Celltype ") + 
  theme(legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0.5, face='bold'),
        axis.text.x = element_text( size = 12, angle = 45, vjust=1, hjust=1, color='black'),
        axis.text.y = element_text( size = 12,  vjust=1, hjust=1, color='black'),
        strip.background = element_blank(),
        strip.text = element_blank()
  ) +

  scale_fill_manual(values=colorset)
p

ggsave(filename, plot = p, width=10 ,height=4.5, 
       dpi=600, limitsize = FALSE)


## trans-eQTL summary
## df_egene dataframe with columns celltype, cis_eGenes, trans_eGenes 
egene_number_plot <- df_egene %>% ggplot(aes(x=cis_eGenes, y=trans_eGenes, 
                        color=celltype)) +
  geom_point(stat = "identity",size=4 ) +
  scale_color_manual(values = colorset) +
  theme_bw() + 
  scale_y_continuous(name = "trans eGenes (n)") +
  scale_x_continuous(name = "cis eGenes (n)", ) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.9, size = 12, color='black'),
        axis.text.y = element_text(angle = 45, size = 12, color='black'),
        axis.title = element_text(size = 12, color='black'),
        legend.position = 'none') +
    ggrepel::geom_text_repel(data = df_egene, aes(label=celltype),
                          force=1, 
                          min.segment.length = unit(0, 'lines'), 
                  box.padding=unit(0.5,'lines'), size=4,
                  direction = "y", hjust = 0, vjust = 1, 
                  color='black' ) 
show(egene_number_plot)

filename <- paste0('~/transQTL_eGenes.pdf')
ggsave(filename, plot = egene_number_plot)

## trans-eQTL upset plot
## transqtl_{celltype}_sig: filter by Bornferroni threshold top assoc.
library('UpSetR')
trabsqtl_share <- merge(transqtl_ast_sig %>% dplyr::select(feature) %>% mutate(Ast=1) %>% unique(), 
                        transqtl_ext_sig %>% dplyr::select(feature) %>% mutate(Ext=1) %>% unique(), 
                        all=TRUE, by = 'feature')
trabsqtl_share <- merge(trabsqtl_share, 
                        transqtl_in_sig %>% dplyr::select(feature) %>% mutate(IN=1) %>% unique(), 
                        all=TRUE, by = 'feature')
trabsqtl_share <- merge(trabsqtl_share, 
                        transqtl_mg_sig %>% dplyr::select(feature) %>% mutate(MG=1) %>% unique(), 
                        all=TRUE, by = 'feature')
trabsqtl_share <- merge(trabsqtl_share, 
                        transqtl_od_sig %>% dplyr::select(feature) %>% mutate(OD=1) %>% unique(), 
                        all=TRUE, by = 'feature')
trabsqtl_share <- merge(trabsqtl_share, 
                        transqtl_opc_sig %>% dplyr::select(feature) %>% mutate(OPC=1) %>% unique(), 
                        all=TRUE, by = 'feature')

trabsqtl_share <- trabsqtl_share %>% column_to_rownames('feature')

trabsqtl_share[is.na(trabsqtl_share)] <- 0

pdf("~/transQTL_eGenes_upset.pdf")
upset(trabsqtl_share, 
      matrix.color="black", nsets = 10, 
      nintersects = 21, point.size=4,
      sets.bar.color=c( "#e31a1c","#ff7f00","#33a02c",
                       "#6a3d9a","#fdbf6f","#1f78b4"
                       )
      )
dev.off()

## trans-eQTL circlo plot
variant_id <- 'rs123'
transqtl_sig_circle_plot <- transqtl_sig[transqtl_sig$variant_id==variant_id,]
transqtl_sig_circle_plot
cis_genename<- 'ABC'

plot1_data  <- data.frame(feature = c('ENSG00000123456', transqtl_sig_circle_plot$feature)) %>% unique()
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= plot1_data$feature, 
                              keytype = "GENEID", columns = c("SYMBOL","GENEID"))
colnames(geneIDs1) = c('name', 'feature')  

plot1_data <- merge(plot1_data, geneIDs1, by ='feature')
plot1_data <- merge(plot1_data, phen_mdata_plot, by ='feature') # phen_mdata_plot includes chr,start,end site for each gene
plot1_data <- data.frame(
  chr = plot1_data$chr,
  start = plot1_data$start,
  end = plot1_data$end,
  name = plot1_data$name
  )
plot1_data

### Define the links for the first plot
plot1_links <- data.frame(
  chr1 = plot1_data[plot1_data$name==cis_genename,c('chr')] %>% as.character(),
  start1 = plot1_data[plot1_data$name==cis_genename,c('start')] %>% as.numeric(), # Approximate center of cis
  end1 = plot1_data[plot1_data$name==cis_genename,c('end')] %>% as.numeric(),
  chr2 = plot1_data[plot1_data$name!=cis_genename,c('chr')] %>% as.character(),
  start2 = plot1_data[plot1_data$name!=cis_genename,c('start')] %>% as.numeric(), # Approximate centers of trans
  end2 = plot1_data[plot1_data$name!=cis_genename,c('end')] %>% as.numeric()
)
### Create the first plot
filename <- paste0('~/transeqtl_circleplot_',cis_genename,'.pdf')
pdf(filename, width = 4, height = 4)
### Create the first plot (page 1 of the PDF)
par(mar = c(2, 2, 2, 2))
circos.par(gap.degree = 1.5, start.degree = 90)
circos.initializeWithIdeogram(
  species = "hg38", 
  plotType = c("ideogram", "labels"),
  chromosome.index = paste0("chr", 1:22), 
)
### Add the connection_height argument to move labels out
circos.genomicLabels(
  plot1_data, 
  cex=0.8,
  labels.column = 4, 
  side = "inside",
  connection_height = mm_h(2),
)
circos.genomicLink(plot1_links[,1:3], plot1_links[,4:6], col = "#000000")
title(paste0("celltype\n",variant_id))
circos.clear()

dev.off()
