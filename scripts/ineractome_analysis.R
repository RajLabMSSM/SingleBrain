suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(LDlinkR))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(valr))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ggrepel))

### Set up the Lifeover function
chain_hg19_hg38 <- import.chain("~/hg19ToHg38.over.chain")
chain_hg38_hg19 <- import.chain("~/hg38ToHg19.over.chain")
liftOverCoord <- function(coord_string, from = "hg19", to = "hg38"){
  if( from == "hg19" & to == "hg38" ){
    chain <- chain_hg19_hg38
  }else{
    stop("only hg19 -> hg38 supported currently")
  }
  coord <- splitCoords(coord_string)
  # lift over assumes chr1 format 
  if( !grepl("chr", coord$chr) ){
    coord$chr <- paste0("chr", coord$chr)
  } 
  # make genomicRanges object
  # liftOver using 
  coord_gr <-  GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(coord$chr),
    ranges = IRanges::IRanges(start=as.numeric(coord$start), end=as.numeric(coord$end)),
    strand = S4Vectors::Rle(rep('*',nrow(coord)))
  )
  # lift over using rtracklayer
  # watch out for duplicate entries
  message(" * lifting over....")
  lifted_over <- rtracklayer::liftOver(coord_gr, chain )
  message(" * lifted!" )
  #return(lifted_over)
  
  stopifnot(length(lifted_over) == length(coord_gr)  )
  coord$chr <- seqnames(lifted_over)
  coord$start <- start(lifted_over)
  coord$end <- end(lifted_over)
  coord$chr <- gsub("chr", "", coord$chr)
  
  lifted_string <- joinCoords(coord, flank = 0)
  return(lifted_string)
}


full_path_gwas <- paste0('~/', gwas_gene,'/',gwas,'_',rsid, '.tsv')
full_path_mm_qtl <- paste0('~/',gwas_gene,'/',qtl,'_',rsid, '.tsv')

### Load the credible SNPs
df_cre <- read_tsv(paste0(path_dir, 'all_COLOC_credible_results_merged_H4_0.5.tsv.gz')) %>% 
  dplyr::filter(GWAS==gwas) %>% 
  dplyr::filter(grepl(gwas_gene, locus)) %>%
  dplyr::filter(QTL==qtl) %>%
  dplyr::filter(grepl(geneid, feature)) %>%
  dplyr::filter(CS >=0.95)

credible_list <- df_cre$snp

gwas_ad <- read_tsv(full_path_gwas) 
gwas_ad$log10p <- -log10(gwas_ad$pvalues)


# Locus zoom plost GWAS and QTL
plot_variants_in_ld_fine <- function(ldlink_results, start, end, rsid, fine) {
  
  ldlink_results <- ldlink_results%>%arrange(PVALUE)
  df_top <- ldlink_results[which(ldlink_results$snp %in% rsid) ,]
  df_fine <- ldlink_results[which(ldlink_results$snp %in% fine) ,]
  df_back <- ldlink_results[!(ldlink_results$snp %in% fine) ,]
  
  
  ldlink_sp <- ggplot(df_back, aes(x=BP, y=PVALUE)) +
    geom_point(aes(colour = R2), binaxis='y', stackdir='center', size = sizes, alpha=0.7) +
    scale_colour_gradientn(colours = c("#1f78b4","#ffff99","#e31a1c"),
                           # scale_colour_distiller(palette = "Spectral",
                           limits = c(0,1), 
                           breaks = c(0,1), 
                           labels = c(0,1) ) +
    geom_point(data=df_fine, aes(x=BP, y=PVALUE, fill='#6a3d9a'), 
               color='#6a3d9a', size = sizes+2, shape=18) +
    # scale_fill_gradientn(colours = c("#1f78b4","#ffff99","#e31a1c"), limits = c(0,1)) +
    xlab(element_blank()) + ylab('-log10(p-value)') +
    scale_x_continuous(breaks = seq(start, end, by = flank)) +
    theme_classic() +
    guides(colour = guide_colourbar(barwidth = 0.5, barheight = 3,
                                    label.position = "right",
                                    ticks  = FALSE, raster =TRUE ),
           fill = "none"
    ) +
    theme(plot.title = element_text(color='black', size = 12, face = 'bold', vjust=0.5),
          axis.title = element_text(color='black',size = 12, vjust=0.5),
          axis.text = element_text(color='black',size = 12,  vjust=0.5),
          legend.text = element_text(color='black',size = 12),
          legend.ticks = element_blank(),
    )
  
  
  
ldlink_sp <- ldlink_sp +
    geom_text_repel(data = df_top, aes(label=snp),
                    force=1,  box.padding=unit(1,'lines'), 
                    direction = "x", hjust = 0.5 ) 
  df_tmp <- df_fine %>% dplyr::filter(PVALUE>4)%>% 
    dplyr::filter(snp!=rsid)
  qtl_lead_snp <- ldlink_results[which.max(ldlink_results$PVALUE), ]
  if(rsid!=qtl_lead_snp$snp){
    # df_qtl <- ldlink_results[!(ldlink_results$snp %in% qtl_rsid) ,]
    ldlink_sp <- ldlink_sp +
      geom_text_repel(data = qtl_lead_snp, aes(label=snp),
                      force=1,  box.padding=unit(1,'lines'), 
                      # direction = "x", 
                      hjust = 0.5 ) 
  }
  if(nrow(df_tmp)>0){
    ldlink_sp <- ldlink_sp +
      geom_text_repel(data = df_tmp %>% dplyr::filter(snp!=qtl_lead_snp$snp), aes(label=snp),
                      force=1,  box.padding=unit(1,'lines'), 
                      # direction = "x", 
                      hjust = 0.5 ) 
  }
  
  
  return(ldlink_sp)
}

### Download LD with LDplinkR
gene_variants_in_ld <- LDproxy(rsid, pop = "EUR", r2d = "r2", token = '')
gene_variants_in_ld <- gene_variants_in_ld[, c(1,7)]
names(gene_variants_in_ld)[names(gene_variants_in_ld) == 'RS_Number'] = 'snp'

### plot the GWAS fine-mapping bar plot
gwas_fine <- fread('~/GWAS_finemapping.csv.gz'),
)
gwas_fine_info <- gwas_fine %>% 
  dplyr::filter(grepl(gwas_gene,Locus)) %>%
  dplyr::filter(FINEMAP.PP>0.1)

gwas_fine_plot <- gwas_fine_info %>% ggplot(aes(x=POS, y=FINEMAP.PP, color=Locus)) +
  geom_point(fill='#6a3d9a' ,size = sizes, shape=18) +
  theme_classic() +
  scale_color_manual(values = c('#6a3d9a','#6a3d9a'  )) +
  xlab(element_blank()) + ylab('PIP') +
  scale_x_continuous(breaks = seq(start_hg19, end_hg19, by = flank)) +
  coord_cartesian(xlim = c(start_hg19, end_hg19), c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  guides(fill = "none") +
  theme(plot.title = element_text(color='black', size = 12, face = 'bold', vjust=0.5),
        axis.title = element_text(color='black',size = 12, vjust=0.5),
        axis.text.y = element_text(color='black',size = 12,  vjust=0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
  ) + geom_hline(yintercept=c(0.95), linetype="solid",
                 color = "gray", size=0.5, )
gwas_fine_plot

### GWAS locuszoom plot
gwas_fine <- gwas_fine %>% #arrange(desc(PIP)) %>%
  dplyr::filter(grepl(gwas_gene,Locus)) %>%
  dplyr::filter((FINEMAP.PP>0.95)) #%>% 

gene_variants_to_plot <- inner_join(gene_variants_in_ld, gwas_ad)
names(gene_variants_to_plot)[names(gene_variants_to_plot) == 'pos'] = 'BP'
names(gene_variants_to_plot)[names(gene_variants_to_plot) == 'log10p'] = 'PVALUE'

if (max(gene_variants_to_plot$PVALUE) > 10 ){
  y_max <- (max(gene_variants_to_plot$PVALUE)+2)
}else{
  y_max <- 12
}

gwas_finemap_list <- gwas_fine$SNP %>% unique()
gwas_finemap_list <-gwas_finemap_list[which(gwas_finemap_list %in% gene_variants_to_plot$snp)]

gwas_ad_fine_plot <- plot_variants_in_ld_fine(gene_variants_to_plot, 
                                              start_hg19, end_hg19, rsid, 
                                              gwas_finemap_list) + 
  ggtitle(paste0(disease,' GWAS : ',gwas_gene)) +
  coord_cartesian(xlim = c(start_hg19, end_hg19), c(0, y_max))
gwas_ad_fine_plot

### Plot the locuszoom with fine-mapping 
qtl_fine <- read_tsv('~/all_results_merged_PIP_0.1.tsv.gz')
qtl_fine <- qtl_fine %>% dplyr::filter((QTL==qtl) & (grepl(geneid, Gene_Ensembl)))  %>%
  dplyr::filter(method==select_method) %>%
  dplyr::filter(PIP >=0.95 )

qtl_variants_to_plot <- inner_join(gene_variants_in_ld, qtl_mm)
names(qtl_variants_to_plot)[names(qtl_variants_to_plot) == 'pos'] = 'BP'
names(qtl_variants_to_plot)[names(qtl_variants_to_plot) == 'log10p'] = 'PVALUE'

qtl_variants_to_plot <- qtl_variants_to_plot[!is.na(qtl_variants_to_plot$PVALUE),]
if (max(qtl_variants_to_plot$PVALUE) > 10 ){
  y_max <- (max(qtl_variants_to_plot$PVALUE)+2)
}else{
  y_max <- 12
}

qtl_lead_snp <- qtl_variants_to_plot[which.max(qtl_variants_to_plot$PVALUE), 'snp']
qtl_finemap_list <-qtl_fine$SNP[which(qtl_fine$SNP %in% qtl_variants_to_plot$snp)]

qtl_mm_fine_plot <- plot_variants_in_ld_fine(qtl_variants_to_plot, 
                                             start, end, 
                                             rsid, qtl_finemap_list) +
  ggtitle(paste0(qtl_name,' eQTL : ',gene_name)) +
  coord_cartesian(xlim = c(start, end), ylim = c(0, y_max)) 
qtl_mm_fine_plot


#FINE BAR
qtl_fine <- read_tsv('~/all_results_merged_PIP_0.1.tsv.gz')
qtl_fine <- qtl_fine %>% dplyr::filter((QTL==qtl) & (grepl(geneid, Gene_Ensembl)))  %>%
  dplyr::filter(method==select_method)
qtl_fine <-qtl_fine[which(qtl_fine$SNP %in% qtl_variants_to_plot$snp),]

qtl_mm_finebar_plot <- qtl_fine %>%  ggplot(aes(x=BP, y=PIP, color=method)) +
  geom_point(fill='#6a3d9a' ,size = sizes, shape=18) +
  theme_classic() +
  scale_color_manual(values = c('#6a3d9a' )) +
  xlab(element_blank()) + ylab('PIP') +
  scale_x_continuous(breaks = seq(start, end, by = flank)) +
  coord_cartesian(xlim = c(start, end), c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  guides(fill = "none") +
  theme(plot.title = element_text(color='black', size = 12, face = 'bold', vjust=0.5),
        axis.title = element_text(color='black',size = 12, vjust=0.5),
        axis.text.y = element_text(color='black',size = 12,  vjust=0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
  ) + geom_hline(yintercept=c(0.95), linetype="solid", 
                 color = "gray", size=0.5, )
qtl_mm_finebar_plot

### Add line for lead GWAS, eQTL, fine-mapping SNP 
snp_loc <- gwas_ad[which(gwas_ad$snp %in% snp_target_list),]
gwas_ad_fine_plot2 <-  gwas_ad_fine_plot +
  geom_vline(xintercept=c( snp_loc$pos), linetype="solid",
             color = "#e31a1c", size=0.1, )
gwas_fine_plot2 <-  gwas_fine_plot +
  geom_vline(xintercept=c(snp_loc$pos), linetype="solid",
             color = "#e31a1c", size=0.1, )
qtl_mm_fine_plot2 <-  qtl_mm_fine_plot +
  geom_vline(xintercept=c(snp_loc_new$end), linetype="solid",
             color = "#e31a1c", size=0.1, )
qtl_mm_finebar_plot2 <-  qtl_mm_finebar_plot +
  geom_vline(xintercept=c(snp_loc_new$end), linetype="solid",
             color = "#e31a1c", size=0.1, )


plot_zoom <- cowplot::plot_grid(gwas_ad_fine_plot2,
                                gwas_fine_plot2,
                                qtl_mm_fine_plot2,
                                qtl_mm_finebar_plot2,
                                align='v', nrow=4, ncol=1,
                                rel_heights = c(4, 1,
                                                4, 1 
                                ))

plot_zoom

ggsave(filename, plot = plot_zoom, units = c("in"),
       width=7 ,height=9, dpi=600,
       limitsize = FALSE)


## Locus zoom crossdot plot 
plot_variants_cross <- function(ldlink_results, gwas_rsid,
                                gwas_fine_list, qtl_fine_list) {
  
  df_gwas_top <- ldlink_results[which(ldlink_results$snp %in% gwas_rsid) ,]
  df_qtl_top <- ldlink_results[which.max(ldlink_results$PVALUE) ,]
  df_gwas_fine <- ldlink_results[which(ldlink_results$snp %in% gwas_fine_list) ,]
  df_qtl_fine <- ldlink_results[which(ldlink_results$snp %in% qtl_fine_list) ,]
  
  df_plot <- ldlink_results[!(ldlink_results$snp %in% gwas_fine_list) ,]
  df_plot <- df_plot[!(df_plot$snp %in% qtl_fine_list) ,]
  
  # ldlink_sp ldlink_results %>% dplyr::filter(snp==snp )%>% unique()
  ldlink_sp <- ggplot(df_plot, aes(x=log10p_gwas, y=PVALUE)) +
    geom_point(aes(colour = R2), binaxis='y', stackdir='center', size = sizes, alpha=0.7) +
    scale_colour_gradientn(colours = c("#1f78b4","#ffff99","#e31a1c"),
                           # scale_colour_distiller(palette = "Spectral",
                           limits = c(0,1), 
                           breaks = c(0,1), 
                           labels = c(0,1) ) +
    xlab('-log10(p-value) GWAS') + ylab('-log10(p-value) eQTL') +
    geom_point(data=df_gwas_fine, aes(x=log10p_gwas, y=PVALUE
    ), fill='#6a3d9a', color='#6a3d9a',
    size = sizes, shape=24) +
    geom_point(data=df_qtl_fine, aes(x=log10p_gwas, y=PVALUE
    ), fill='#6a3d9a', color='#6a3d9a',
    size = sizes, shape=25) +  
    theme_bw() +
    ggtitle(paste0(disease,'-',gwas_gene,'-',file_subname,'-',gene_name)) +
    theme( axis.text.x = element_text(colour = "black", size=12),
           axis.title.x = element_text(colour = "black",size=12),
           axis.text.y = element_text(colour = "black",size=12),
           axis.title.y = element_text(colour = "black",size=12),
           legend.text = element_text(color='black',size = 12),
           legend.ticks = element_blank(),
    )+
    guides(colour = guide_colourbar(barwidth = 0.5, barheight = 3,
                                    label.position = "right",
                                    ticks  = FALSE, raster =TRUE ),
           fill = "none"
    )
  
  df_text <- rbind(df_gwas_top, df_qtl_top)
  df_text <- rbind(df_text, df_gwas_fine)
  df_text <- rbind(df_text, df_qtl_fine) %>% unique()
  
  ldlink_sp <- ldlink_sp +
    geom_text_repel(data = df_text, aes(label=snp),
                    color='black',
                    size=4,
                    force=1,  
                    box.padding=unit(1,'lines'),
                    # direction = "x", 
                    hjust = 0.5 )
  
  
  return(ldlink_sp)
  
}

cross_gwas <- gene_variants_to_plot %>% dplyr::select(snp,  log10p_gwas=PVALUE)
all_res_only_tmp <- inner_join(cross_gwas, qtl_variants_to_plot)

plot_cross <- plot_variants_cross(all_res_only_tmp, rsid,
                                  gwas_finemap_list, qtl_finemap_list)
plot_cross

ggsave(filename, plot = plot_cross, units=c('in'),
       width = 6.5, height = 6, dpi = 600,
       limitsize = FALSE
)



## Interactome analysis

### SNP and Nott et al. dot plot
snp_loc <- gene_variants_to_plot %>% dplyr::filter(pvalues <= 1e-4)
finemap_list <- c(qtl_finemap_list, gwas_finemap_list,  rsid, qtl_lead_snp, credible_list) #rsid_tmp,
snp_loc <- snp_loc[which(snp_loc$snp %in% finemap_list),]

df_snp <- snp_loc[,c('snp','R2')]

#For hg38 GWAS
snp_loc$chr <- paste0("chr", snp_loc$chr)
grange_loc <-
    GenomicRanges::GRanges(seqnames = snp_loc$chr,
            ranges = IRanges::IRanges(start = snp_loc$BP - 1, end = snp_loc$BP),
            SNP = snp_loc$snp)#, GWAS = snp_loc$GWAS, label = snp_loc$label )

message(" * lifting over....")
lifted_over <- rtracklayer::liftOver(grange_loc, chain_hg38_hg19 )
message(" * lifted!" )
snp_loc_rsid <- GenomeInfoDb::as.data.frame(lifted_over)

names(snp_loc_rsid)[names(snp_loc_rsid) == 'SNP'] = 'snp'
names(snp_loc_rsid)[names(snp_loc_rsid) == 'seqnames'] = 'chr'

#For hg19 GWAS
snp_loc_rsid <- snp_loc
snp_loc_rsid$chr <- paste0("chr", snp_loc_rsid$chr)
snp_loc_rsid$end <- snp_loc_rsid$BP
snp_loc_rsid$start <- snp_loc_rsid$end - 1

snp_loc_rsid <- merge(snp_loc_rsid, df_snp)

df_snp_gwas_lead <- snp_loc_rsid[which(snp_loc_rsid$snp %in% rsid),]
df_snp_gwas_lead$label <- 'GWAS Lead SNP'

df_snp_gwas_fine <- snp_loc_rsid[which(snp_loc_rsid$snp %in% gwas_finemap_list),]
df_snp_gwas_fine$label <- 'GWAS Fine-mapping'

df_snp_qtl_lead <- snp_loc_rsid[which(snp_loc_rsid$snp %in% qtl_lead_snp),] 
df_snp_qtl_lead$label <- 'eQTL Lead SNP'

df_snp_qtl_fine <- snp_loc_rsid[which(snp_loc_rsid$snp %in% qtl_finemap_list),] 
df_snp_qtl_fine$label <- 'eQTL Fine-mapping'

df_snp_credible <- snp_loc_rsid[which(snp_loc_rsid$snp %in% credible_list),]
df_snp_credible$label <- 'Credible SNPs'

locus_res <- rbind(df_snp_gwas_lead, df_snp_gwas_fine)
locus_res <- rbind(locus_res, df_snp_qtl_lead)
locus_res <- rbind(locus_res, df_snp_qtl_fine)
locus_res <- rbind(locus_res, df_snp_credible)

locus_res$R2 <- signif(locus_res$R2, digits = 2)

locus_res <- locus_res %>% arrange(end)
snp_seq <- locus_res$snp %>% unique()

### Plot the SNP by category
r2_plot <- locus_res %>%
  mutate(snp = factor(snp, levels = snp_seq)) %>%
  mutate(label = factor(label, levels = c('GWAS Lead SNP', 'GWAS Fine-mapping', 
                                          'eQTL Lead SNP', 'eQTL Fine-mapping',
                                          'Credible SNPs'
  ))) %>%
  mutate(label = forcats::fct_rev(label)) %>%
  ggplot(aes(x = snp, y = label ) ) + 
  geom_point(aes(colour = R2), size = 10, position = position_nudge(x = 0)) +
  geom_text(aes(label = R2), nudge_x = 0, size = 3) +
  theme_bw() +
  scale_x_discrete(position = 'top', guide = guide_axis(n.dodge = 2)) +
  scale_colour_gradientn(colours = c("#1f78b4","#ffff99","#e31a1c"), 
                         limits = c(0,1), 
                         breaks = c(0,1), 
                         labels = c(0,1) ) +
  labs(y = "", x = "",
       title = paste("Locus:", gwas_gene, "Gene: ", gene_name), 
  ) +
  guides(colour = guide_colourbar(barwidth = 0.5, barheight = 3,
                                  label.position = "right",
                                  ticks  = FALSE, raster =TRUE )) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black", size=10), 
        axis.text.y = element_text(colour = "black", size=12),
        axis.ticks = element_line(colour = "black"),
        legend.text = element_text(color='black',size = 12),
        legend.ticks = element_blank(),
  )

r2_plot

### Check the SNP with Nott et al. 2019 enhancer and promoter sites
nott_overlap <- function(locus_res, set, return_coordinates = FALSE){
  snp_df <- dplyr::filter(locus_res, !is.na(chr), !is.na(start)) 
  snp_gr <- GenomicRanges::GRanges(seqnames = snp_df$chr, ranges = IRanges::IRanges(start = snp_df$start, end = snp_df$end), snp = snp_df$snp)
  
  df <- nott_2019_interactome(set)
  df_gr <- GenomicRanges::GRanges(seqnames = df$chr, ranges = IRanges::IRanges(start = df$start, end = df$end))
  overlap <- GenomicRanges::findOverlaps(snp_gr, df_gr)
  set_string <- gsub(" ", "_", set)
  locus_res[[set_string]] <- locus_res$snp %in% snp_gr[overlap@from]$snp
  
  if( return_coordinates == TRUE){
    set_df <- as.data.frame(df_gr[overlap@to])
    coords <- paste0(gsub("chr", "", set_df$seqnames), ":", set_df$start, "-", set_df$end )
    coord_name <- paste0(set_string, "_coord_hg19")
    
    locus_res[[coord_name]] <- NA
    locus_res[[coord_name]][overlap@from] <- coords
  }
  return(locus_res)
}

add_nott_overlap <- function(locus_res, return_coords = FALSE){
  if(is.null(locus_res)){return(NULL)}
  locus_res %>%
    nott_overlap("Microglia enhancers", return_coords) %>%
    nott_overlap("Microglia promoters", return_coords) %>%
    nott_overlap("Neuronal enhancers", return_coords) %>%
    nott_overlap("Neuronal promoters", return_coords) %>%
    nott_overlap("Oligo enhancers", return_coords) %>%
    nott_overlap("Oligo promoters", return_coords) %>%
    nott_overlap("Astrocyte enhancers", return_coords) %>%
    nott_overlap("Astrocyte promoters", return_coords)
}

locus_res <- add_nott_overlap(locus_res)
locus_res


to_plot <- locus_res %>%
  dplyr::select(snp, start, ends_with("enhancers"), ends_with("promoters")) %>%
  distinct() %>%
  tidyr::gather( key = "feature",
                 value = "overlap", 
                 -snp, -start) %>%
  tidyr::separate(feature, into = c("cell_type", "seq_type"), sep = "_") %>%
  dplyr::mutate(cell_type = case_when(cell_type == "Neuronal" ~"Neurons",
                                      cell_type == "Oligo" ~ "Oligos",
                                      cell_type == "Astrocyte" ~ "Astrocytes",
                                      TRUE ~ cell_type)) %>%
  mutate(cell_type = factor(cell_type, levels = rev(c("Astrocytes", "Microglia", "Neurons", "Oligos"))))

nott_plot <- to_plot %>% 
  mutate(snp = factor(snp, levels = unique(locus_res$snp))) %>%
  ggplot(aes(x = snp, y = cell_type)) + 
  geom_point(aes(colour = seq_type,alpha = overlap), size = 10 ) +
  labs(y = "", x = "", colour = "Overlapping") + 
  theme_bw() +
  scale_alpha_manual(values = c(0,1)) +
  scale_colour_manual(values = c("enhancers" = "goldenrod", "promoters" = "darkviolet")) +
  theme(panel.grid = element_blank() ) +
  guides(alpha = FALSE) +
  theme(axis.text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black", size=12),
        axis.text.x = element_text(colour = "black", size=12, angle = 45, 
                                   vjust = 1, hjust = 1),
        axis.ticks = element_line(colour = "black"),
        legend.text = element_text(color='black',size = 12),
        legend.ticks = element_blank(),
  )
nott_plot

### Merge R2 and Nott et al. plot
overlap_result <- cowplot::plot_grid(r2_plot, 
                                     nott_plot,
                                     align='v', nrow=2, ncol=1,
                                     rel_heights = c(7,  
                                                     6))
overlap_result

ggsave(filename, plot = overlap_result, units = c("in"),
       width=8 ,height=6, dpi=600,
       limitsize = FALSE)



# Interactome analysis Part2
# Functions
## SNP plot
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 111)

db.gr <- ensembldb::transcripts(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86) %>%
  data.table::as.data.table() %>%
  dplyr::filter( tx_biotype == "protein_coding") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::slice_max(width, n = 1) %>%
  distinct(gene_id, .keep_all = TRUE)

snp_plot <- function(df_snp_loc, start, end, snp_width = 50){
  require(ggplot2)
  require(GenomicRanges)
  require(ggbio)
  message(length(df_snp_loc$snp))
  grange_loc <-  
    GenomicRanges::GRanges(seqnames = df_snp_loc$chr, 
                           ranges = IRanges::IRanges(start = df_snp_loc$end - snp_width, end = df_snp_loc$end), 
                           SNP = df_snp_loc$snp,
                           label = df_snp_loc$chr)
  grange_loc$SNP <- factor(grange_loc$SNP, levels = df_snp_loc$snp)
  message(length(grange_loc$SNP))

  plot <- 
    ggplot() + 
    ggbio::geom_rect(data = grange_loc, aes(group = label)  ) +
    theme_classic() + 
    coord_cartesian(xlim = c(start, end)) +
    theme(axis.text = element_blank(),
      axis.ticks.x = element_blank(),
    )  

  
  return(plot)
}

## NOTT Plot
nott_epigenome_plot <- function(df_snp_loc, start, end, color){
  plot <- ggplot() + geom_rect(df_snp_loc, mapping=aes(xmin=start,xmax=end, ymin=0, ymax=1, fill = label)) +
    scale_fill_manual(values = c(color)) +
    coord_cartesian(xlim = c(start, end)) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"), 
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.title=element_blank())
  return(plot)
}

nott_plac_plot <- function(plac, snp_loc, start_decoy, end_decoy, color){
  
  snp_loc_start <- start_decoy
  snp_loc_end <- end_decoy
  
  plac_loc <- dplyr::filter(plac_df, start1 >= snp_loc_start & end2 <= snp_loc_end )
  
  grange_loc <-  
    GenomicRanges::GRanges(seqnames = snp_loc$chr, 
                           ranges = IRanges::IRanges(start = snp_loc$end - 1, end = snp_loc$end), 
                           SNP = snp_loc$snp)
  
  plac_start_gr <- GenomicRanges::GRanges( seqnames = plac_loc$chr1,
                                           ranges = IRanges::IRanges(start = plac_loc$start1, end = plac_loc$end1) )
  
  plac_end_gr <- GenomicRanges::GRanges( seqnames = plac_loc$chr2,
                                         ranges = IRanges::IRanges(start = plac_loc$start2, end = plac_loc$end2) )
  # find overlaps
  end_overlaps <- GenomicRanges::findOverlaps(grange_loc, plac_end_gr)
  start_overlaps <- GenomicRanges::findOverlaps(grange_loc, plac_start_gr)
  all_overlaps <- unique(c( S4Vectors::subjectHits(end_overlaps), S4Vectors::subjectHits(start_overlaps)))
  
  plac_loc$overlap <- FALSE
  plac_loc$overlap[all_overlaps] <- TRUE
  
  plac_loc <- plac_loc %>% arrange(overlap)
  ## make plot  
  ggplot() + 
    ggbio::geom_arch(data = plac_loc, aes(x = end1, xend = start2, colour = overlap, alpha = overlap) ) + theme_classic() + labs(x = "", y = "") + 
    ggbio::geom_rect(data = plac_loc, aes(xmin = start1, xmax = end1, fill = overlap, alpha = overlap ), ymin = 0, ymax = 1, colour = NA ) +
    ggbio::geom_rect(data = plac_loc, aes(xmin = start2, xmax = end2, fill = overlap, alpha = overlap), ymin = 0, ymax = 1, colour = NA) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank() ,
          axis.ticks.x = element_blank() ,
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          legend.title=element_blank()
    ) +
    scale_fill_manual(values = c("darkgray", color)) +
    scale_colour_manual(values = c("gray", color)) +
    scale_alpha_manual(values = c(0.1, 1)) +
    guides(colour = FALSE, alpha = FALSE) 
  
  
  
}

## Gene Plot
gene_plot <- function(chr_decoy, start_decoy, end_decoy){
  
  snp_loc_start <- start_decoy
  snp_loc_end <- end_decoy
  snp_gr <- GenomicRanges::GRanges(seqnames = paste0('chr', chr_decoy, sep=''), 
                                   ranges = IRanges::IRanges(start = snp_loc_start, end = snp_loc_end))
  GenomeInfoDb::seqlevelsStyle(snp_gr) <- "NCBI"
  
  db_loc <- subset(db.gr,
                   seqnames == chr_decoy &
                     (start >= snp_loc_start | end <= snp_loc_end) ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(db_loc) <- "NCBI"
  
  edb <-  ensembldb::addFilter( EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                AnnotationFilter::TxIdFilter(db_loc$tx_id) )
  
  plot <- ggplot() +  scale_x_continuous(breaks = seq(snp_loc_start, 
                                                      snp_loc_end, 
                                                      by = (snp_loc_end-snp_loc_start)/2)) +
    coord_cartesian(xlim = c(snp_loc_start, snp_loc_end)) + 
    ggbio::geom_alignment(edb, which = snp_gr,
                          names.expr = "gene_name",
                          aes(group = gene_name,
                              fill='gene',
                              arrow.rate = 0.5,
                              length = unit(0.02, "cm"),
                              show.legend=FALSE)) + 
    scale_fill_manual(values = c("#7A7A7A")) +
    theme_classic() + 
    theme(legend.title=element_blank(),
          axis.text = element_text(color='black',size = 12,  vjust=0.5), )
  
  return(plot)
}

snp_loc <- locus_res %>% dplyr::select(snp, chr, end) %>% distinct()

start_hg19 <- min(snp_loc$end) - width
end_hg19 <- max(snp_loc$end) + width + extra_flank
chr <- unique(snp_loc$chr)

grange_loc <-  
  GenomicRanges::GRanges(seqnames = snp_loc$chr, 
                         ranges = IRanges::IRanges(start = snp_loc$end - 1, end = snp_loc$end), 
                         SNP = snp_loc$snp)

# message(" * lifting over....")
lifted_over <- rtracklayer::liftOver(grange_loc, chain_hg19_hg38 )
# message(" * lifted!" )
snp_loc_new <- GenomeInfoDb::as.data.frame(lifted_over)
start <- min(snp_loc_new$end) - width
end <- max(snp_loc_new$end) + width + extra_flank

SNP_plot <- snp_plot(snp_loc, start_hg19, end_hg19)


cell_promoter <- nott_2019_interactome('Microglia promoters')
peaks <- cell_promoter %>% dplyr::filter(chrom == chr, start >= start_hg19)
peaks <- peaks %>% dplyr::filter(chrom == chr, end <= end_hg19)
peaks$label <- 'Microglia_promoter'
microglia_promoter_plot <- nott_epigenome_plot(peaks, start_hg19, end_hg19, "#6a3d9a")
microglia_promoter_plot

cell_enhancer <- nott_2019_interactome('Microglia enhancers')
peaks <- cell_enhancer %>% dplyr::filter(chrom == chr, start >= start_hg19)
peaks <- peaks %>% dplyr::filter(chrom == chr, end <= end_hg19)
peaks$label <- 'Microglia_enhancer'
Microglia_enhancer_plot <- nott_epigenome_plot(peaks, start_hg19, end_hg19, "#6a3d9a")
Microglia_enhancer_plot

snp_loc_plac <- gene_variants_to_plot %>% dplyr::filter(pvalues < 1e-04)
finemap_list <- c(qtl_finemap_list, gwas_finemap_list,  rsid, qtl_lead_snp, credible_list) 
snp_loc_plac <- snp_loc_plac[which(snp_loc_plac$snp %in% finemap_list),]
names(snp_loc_plac)[names(snp_loc_plac) == 'BP'] = 'end'

snp_loc_new <- snp_loc_plac
snp_loc_new$chr <- paste0('chr',snp_loc_new$chr)

plac_df <- read.table('/Users/jin/Desktop/Jin/brain-cell-type-peak-files/PLACseq/Microglia.5k.2.peaks.bedpe',
                      header = T)

Microglia_plac_plot <- nott_plac_plot(plac_df, 
                                      snp_loc_new, 
                                      start_hg19, 
                                      end_hg19,
                                      '#6a3d9a')  + 
  coord_cartesian(xlim = c(start_hg19, end_hg19))
Microglia_plac_plot




gene_transcripts <- gene_plot(gsub('chr','',chr), start, end) 
gene_transcripts

plot_epi <- cowplot::plot_grid( SNP_plot,
                                microglia_promoter_plot,
                                Microglia_enhancer_plot,
                                Microglia_plac_plot,
                                gene_transcripts,
                                align='v', nrow=5, ncol=1,
                                rel_heights = c(1, 
                                                1.5,
                                                1.5, 
                                                1.5,
                                                2
                                ))

plot_epi

ggsave(filename, plot = plot_epi, units = c("in"),
       width=8 ,height=5, dpi=600,
       limitsize = FALSE)




# End


