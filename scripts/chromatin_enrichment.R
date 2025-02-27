
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(rtracklayer))
library(GenomicRanges)

# Enrichment TSS
GENCODE <- rtracklayer::import("~/gencode.v38.primary_assembly.annotation.gtf")

GENCODE_TSS <- GENCODE %>% as.data.frame() %>%
  dplyr::filter(type=='gene')%>%
  dplyr::select(seqnames, start, end, strand, gene_id)

GENCODE_TSS$TSS <- GENCODE_TSS$start
GENCODE_TSS$TSS <- ifelse(grepl("-", GENCODE_TSS$strand), GENCODE_TSS$end, GENCODE_TSS$start)

mmqtl_tss <- read_tsv("eqtl_egene_pos.tsv")
mmqtl_tss <- merge(mmqtl_tss, GENCODE_TSS, by='feature')
mmqtl_tss$distance <- mmqtl_tss$TSS - mmqtl_tss$pos

plot_tss <- mmqtl_tss %>%
  ggplot( aes(x=distance, group=celltype, color=celltype)) +
  geom_density(alpha=1) +
  scale_color_manual(values=colorset2) +
  theme_classic() +
  labs(x='Distance (kbp)', y='Density') +
  scale_x_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3, sep = ""),
                     limits = c(-1500000,1500000)
  )+
  scale_y_continuous() +
  facet_wrap(~celltype, scale="free", ncol = 4) +
  theme(axis.title = element_text(color='black', size=12),
        axis.text.y = element_text(color='black', size=10),
        axis.text.x = element_text(color='black', size=10),
        legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color='black', size=12),
  ) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "grey", size=0.5,)

plot_tss

ggsave(filename, plot = plot_tss, 
       width=12 ,height=5, 
       dpi=600, limitsize = FALSE)


# Enrichiment Enhancer and promoters

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
res

gencode <- rtracklayer::import("~/gencode.v38.primary_assembly.annotation.gtf")
gencode_genes <- gencode[ gencode$type == "gene"]

  gene_meta <- tibble(gene_id = gencode_genes$gene_id,
                      gene_name = gencode_genes$gene_name,
                      chr = as.character(seqnames(gencode_genes)),
                      gene_start = start(gencode_genes),
                      gene_end = end(gencode_genes),
                      strand = as.character(strand(gencode_genes) )
                      )
  
res_gene <- map(res, ~{
    left_join(.x , gene_meta, by = c("feature" = "gene_id", "chr") )
  })

chain_hg19_hg38 <- import.chain(paste0(path_chain, "hg19ToHg38.over.chain"))
chain_hg38_hg19 <- import.chain(paste0(path_chain, "hg38ToHg19.over.chain"))
liftOverCoord <- function(dataset, from = "hg19", to = "hg38"){
  if( from == "hg19" & to == "hg38" ){
    chain <- chain_hg19_hg38
  }else if( from == "hg38" & to == "hg19" ){
    chain <- chain_hg38_hg19
  }else{
    stop("only hg19 -> hg38 supported currently")
  }
  
  # lift over using rtracklayer
  # watch out for duplicate entries
  message(" * lifting over....")
  lifted_over <- rtracklayer::liftOver(dataset, chain )
  message(paste0(" * lifted! ",from,' to ',to))
  #return(lifted_over)
  
  snp_loc_new <- GenomeInfoDb::as.data.frame(lifted_over)
  dataset2 <- GenomicRanges::GRanges(seqnames = snp_loc_new$seqnames, 
                                     ranges = IRanges::IRanges(start = snp_loc_new$start, end = snp_loc_new$end))
  
  return(dataset2)
}

### Nott et al. 2019
ast_promoter <- rtracklayer::import(here::here(paste0(path_ext,'/Astrocyte_promoters.bed' )), 
                                    format = "BED")
ast_promoter2 <- liftOverCoord(ast_promoter, from = "hg19", to = "hg38")
ast_enhancer <- rtracklayer::import(here::here(paste0(path_ext,'/Astrocyte_enhancers.bed' )), 
                                    format = "BED")
ast_enhancer2 <- liftOverCoord(ast_enhancer, from = "hg19", to = "hg38")

mg_promoter <- rtracklayer::import(here::here(paste0(path_ext,'/Microglia_promoters.bed' )), 
                                   format = "BED")
mg_promoter2 <- liftOverCoord(mg_promoter, from = "hg19", to = "hg38")
mg_enhancer <- rtracklayer::import(here::here(paste0(path_ext,'/Microglia_enhancers.bed' )), 
                                   format = "BED")
mg_enhancer2 <- liftOverCoord(mg_enhancer, from = "hg19", to = "hg38")

od_promoter <- rtracklayer::import(here::here(paste0(path_ext,'/Oligo_promoters.bed' )), 
                                   format = "BED")
od_promoter2 <- liftOverCoord(od_promoter, from = "hg19", to = "hg38")
od_enhancer <- rtracklayer::import(here::here(paste0(path_ext,'/Oligo_enhancers.bed' )), 
                                   format = "BED")
od_enhancer2 <- liftOverCoord(od_enhancer, from = "hg19", to = "hg38")

neu_promoter <- rtracklayer::import(here::here(paste0(path_ext,'/Neuronal_promoters.bed' )), 
                                    format = "BED")
neu_promoter2 <- liftOverCoord(neu_promoter, from = "hg19", to = "hg38")
neu_enhancer <- rtracklayer::import(here::here(paste0(path_ext,'/Neuronal_enhancers.bed' )), 
                                    format = "BED")
neu_enhancer2 <- liftOverCoord(neu_enhancer, from = "hg19", to = "hg38")

# calculate distances for everything
dist_df <- map_df(res_gene, ~{
  pos_df <- tibble(gene_name= .x$gene_name. ,gene_id = .x$feature, 
                   variant_id= .x$variant_id,  snp_pos = .x$pos, snp_chr = .x$chr,
                   gene_start = .x$gene_start, gene_end = .x$gene_end, strand = .x$strand, 
                   qval = .x$qval, sig = .x$qval < 0.05) %>%
    mutate(dist =  case_when( 
      strand == "+" ~ snp_pos - gene_start, 
      strand == "-" ~ gene_end - snp_pos
    ) )  %>%
    mutate(within_promoter = strand == "+" & snp_pos >= (gene_start - 10000) & snp_pos < (gene_start + 1000) |
             strand == "-" & snp_pos <= (gene_end + 10000) & snp_pos >= (gene_end - 1000)
    ) %>%
    mutate(within_body = (snp_pos > gene_start) & (snp_pos < gene_end) )
  # make Grange for SNPs
  lead_snp_gr <- GRanges(seqnames = .x$chr, ranges = IRanges(
    start = .x$pos - 1, 
    end = .x$pos), 
    id = .x$variant_id)
  # calculate distance to nearest ATAC peak and exon
  pos_df$ast_enh <- as.data.frame(distanceToNearest(lead_snp_gr, subject = ast_enhancer2))$distance
  pos_df$ast_pro <- as.data.frame(distanceToNearest(lead_snp_gr, subject = ast_promoter2))$distance
  pos_df$mg_enh <- as.data.frame(distanceToNearest(lead_snp_gr, subject = mg_enhancer2))$distance
  pos_df$mg_pro <- as.data.frame(distanceToNearest(lead_snp_gr, subject = mg_promoter2))$distance
  pos_df$od_enh <- as.data.frame(distanceToNearest(lead_snp_gr, subject = od_enhancer2))$distance
  pos_df$od_pro <- as.data.frame(distanceToNearest(lead_snp_gr, subject = od_promoter2))$distance
  pos_df$neu_enh <- as.data.frame(distanceToNearest(lead_snp_gr, subject = neu_enhancer2))$distance
  pos_df$neu_pro <- as.data.frame(distanceToNearest(lead_snp_gr, subject = neu_promoter2))$distance
  return(pos_df)
}, .id = "dataset") %>%
  tidyr::separate(col = dataset, into = c("reference", "set"), sep = "_", extra = "merge", remove = FALSE) %>%
  drop_na() # genes that couldn't be matched have NA for start and end

dist_df

qtl_level <- qtl_summary %>% arrange(desc(perc_features))
phenotype <- qtl_level$reference

## plot odds ratios for each QTL set
fisher_test_plot <- function(res, title = "", xlim = NULL){
  p <- res %>%
    mutate(plabel = ifelse( padj < 0.05, "*", "")) %>%
    # mutate(reference = forcats::fct_rev(reference)) %>%
    ggplot( aes(x = reference, y = estimate, colour = reference )) + 
    geom_errorbar(aes(ymin = conf.low, 
                      ymax = conf.high), 
                  width = 0.25, color='black',
                  position = position_dodge(width = 1), 
                  show.legend = FALSE ) + 
    geom_point(size=3 ,
               position = position_dodge(width = 1) ) + 
    facet_grid(. ~external, space='free', scales = 'free', 
               switch = 'x') + 
    scale_color_manual(values = colorset2) +
    # coord_flip() +
    geom_hline(yintercept =  1, linetype = 3) +
    geom_text(aes(x = reference, y = 0.5, label = plabel), 
              position = position_dodge(width = 0.5), 
              size = 8, color='black',
              show.legend = FALSE) +
    labs(y = "Odds ratio", title = title, x = "") + 
    theme_jh() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x= element_blank(),
          # panel.spacing.y = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(colour = "black", size=12),
          strip.placement = "bottom",
          legend.position = "none",
          legend.box.margin = margin(c(0,0,0,0)),
          legend.margin = margin(c(0,0,0,0)),
          #strip.placement = "outside",
          panel.spacing.y = unit(x = 0, units = "points"), 
          plot.title = element_text(colour = "black", 
                                    size=12, face='bold'),
    )
  
  if(! is.null(xlim)){p <- p + scale_y_continuous(limits = xlim) }
  return(p)
}

theme_bj <- function () { 
  theme_bw(base_size=12, base_family="Helvetica") %+replace% 
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(),
      axis.ticks = element_line(colour = "black"),
      #text = element_text(color = "black"), 
      strip.text = element_text(color = "black", size=12),
      axis.text = element_text(colour = "black", size=12),
      plot.text = element_text(colour = "black", size=12),
      panel.background  = element_blank(),
      plot.background = element_rect(fill="white", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA), 
      legend.text = element_text(size = 7)
    )
}

## Enhancer
ast_enh_res2 <- dist_df %>% #dplyr::filter(sig==TRUE) %>%
  split(.$dataset) %>%
  map_df( ~{
    
  }, .id = "set") %>% 
  separate(col = "set", into = c("reference", "phenotype"), sep = "_", extra = "merge", remove = FALSE) 

ast_enh_res2$reference <- factor(ast_enh_res2$reference, levels = phenotype)
ast_enh_res2

mg_enh_res2 <- dist_df %>% #dplyr::filter(sig==TRUE) %>%
  split(.$dataset) %>%
  map_df( ~{
    broom::tidy(fisher.test(table( .x$mg_enh < 1 & .x$mg_enh > -1, sig = .x$sig)))
  }, .id = "set") %>% 
  separate(col = "set", into = c("reference", "phenotype"), sep = "_", extra = "merge", remove = FALSE) 

mg_enh_res2$reference <- factor(mg_enh_res2$reference, levels = phenotype)
mg_enh_res2

od_enh_res2 <- dist_df %>% #dplyr::filter(sig==TRUE) %>%
  split(.$dataset) %>%
  map_df( ~{
    broom::tidy(fisher.test(table( .x$od_enh < 1 & .x$od_enh > -1, sig = .x$sig)))
  }, .id = "set") %>% 
  separate(col = "set", into = c("reference", "phenotype"), sep = "_", extra = "merge", remove = FALSE) 

od_enh_res2$reference <- factor(od_enh_res2$reference, levels = phenotype)
od_enh_res2

neu_enh_res2 <- dist_df %>% #dplyr::filter(sig==TRUE) %>%
  split(.$dataset) %>%
  map_df( ~{
    broom::tidy(fisher.test(table( .x$neu_enh < 1 & .x$neu_enh > -1, sig = .x$sig)))
  }, .id = "set") %>% 
  separate(col = "set", into = c("reference", "phenotype"), sep = "_", extra = "merge", remove = FALSE) 

neu_enh_res2$reference <- factor(neu_enh_res2$reference, levels = phenotype)
neu_enh_res2

ast_enh_res2$external <- 'Astrocytes'
mg_enh_res2$external <- 'Microglia'
od_enh_res2$external <- 'Oligodendrocytes'
neu_enh_res2$external <- 'Neurons'


nott_enh_res <- rbind(ast_enh_res2, mg_enh_res2,
                      od_enh_res2, neu_enh_res2
)

nott_enh_res <- nott_enh_res %>% group_by(external) %>%
  mutate(padj = p.adjust(p.value, method = "bonferroni"))

nott_enh_res$external <- factor(nott_enh_res$external, levels = c('Microglia','Astrocytes', 
                                                                  'Oligodendrocytes','Neurons'
))
# nott_res

nott_enh_plot2 <- fisher_test_plot(nott_enh_res, title = "eQTL lead SNPs enrichements of enhancer peaks")
nott_enh_plot2

ggsave(filename, plot = nott_enh_plot2, 
       width=8 ,height=4, 
       dpi=600, limitsize = FALSE)

ast_enh_res2 <- dist_df %>% #dplyr::filter(sig==TRUE) %>%
  split(.$dataset) %>%
  map_df( ~{
    
  }, .id = "set") %>% 
  separate(col = "set", into = c("reference", "phenotype"), sep = "_", extra = "merge", remove = FALSE) 

ast_enh_res2$reference <- factor(ast_enh_res2$reference, levels = phenotype)
ast_enh_res2


## Promoter
ast_pro_res2 <- dist_df %>%
  split(.$dataset) %>%
  map_df( ~{
    broom::tidy(fisher.test(table( .x$ast_pro < 1 & .x$ast_pro > -1, sig = .x$sig)))
  }, .id = "set") %>% 
  separate(col = "set", into = c("reference", "phenotype"), sep = "_", extra = "merge", remove = FALSE) 

ast_pro_res2$reference <- factor(ast_pro_res2$reference, levels = phenotype)
ast_pro_res2

mg_pro_res2 <- dist_df %>%
  split(.$dataset) %>%
  map_df( ~{
    broom::tidy(fisher.test(table( .x$mg_pro < 1 & .x$mg_pro > -1, sig = .x$sig)))
  }, .id = "set") %>% 
  separate(col = "set", into = c("reference", "phenotype"), sep = "_", extra = "merge", remove = FALSE) 

mg_pro_res2$reference <- factor(mg_pro_res2$reference, levels = phenotype)
mg_pro_res2

od_pro_res2 <- dist_df %>%
  split(.$dataset) %>%
  map_df( ~{
    broom::tidy(fisher.test(table( .x$od_pro < 1 & .x$od_pro > -1, sig = .x$sig)))
  }, .id = "set") %>% 
  separate(col = "set", into = c("reference", "phenotype"), sep = "_", extra = "merge", remove = FALSE) 

od_pro_res2$reference <- factor(od_pro_res2$reference, levels = phenotype)
od_pro_res2

neu_pro_res2 <- dist_df %>%
  split(.$dataset) %>%
  map_df( ~{
    broom::tidy(fisher.test(table( .x$neu_pro < 1 & .x$neu_pro > -1, sig = .x$sig)))
  }, .id = "set") %>% 
  separate(col = "set", into = c("reference", "phenotype"), sep = "_", extra = "merge", remove = FALSE) 

neu_pro_res2$reference <- factor(neu_pro_res2$reference, levels = phenotype)
neu_pro_res2

ast_pro_res2$external <- 'Astrocytes'
mg_pro_res2$external <- 'Microglia'
od_pro_res2$external <- 'Oligodendrocytes'
neu_pro_res2$external <- 'Neurons'

nott_pro_res <- rbind(ast_pro_res2, mg_pro_res2,
                      od_pro_res2, neu_pro_res2
)

nott_pro_res <- nott_pro_res %>% #group_by(external) %>%
  mutate(padj = p.adjust(p.value, method = "bonferroni"))

nott_pro_res$external <- factor(nott_pro_res$external, 
                                levels = c('Microglia', 'Astrocytes',
                                           'Oligodendrocytes','Neurons'
))
# nott_res
nott_pro_plot2 <- fisher_test_plot(nott_pro_res, title = "eQTL lead SNPs enrichements of promoter peaks")
nott_pro_plot2

ggsave(filename, plot = nott_pro_plot2, 
       width=8 ,height=4, 
       dpi=600, limitsize = FALSE)





# End




