suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))


all_res <- read_tsv('~/MPRA_consensus_table.tsv')
sig_res <- all_res %>% dplyr::filter(significance==TRUE)

snp_target_list <- c(qtl_finemap_list, gwas_finemap_list,  rsid, qtl_lead_snp, credible_list)
snp_target_list <- snp_target_list[snp_target_list %in%sig_res$variant_id] %>% unique()

plt_res <- all_res
plt_res$log10p <- -log10(plt_res$pval)
plt_res$sig <- 'No'
plt_res[plt_res$variant_id %in% sig_res$variant_id, 'sig'] = 'Yes'
plt_res[plt_res$variant_id %in% snp_target_list ,'sig'] = 'CredibleSet'
plt_res$sig <- factor(plt_res$sig, levels = c('CredibleSet', 'Yes','No'))

plt_volcano <- ggplot(plt_res, aes(x=logFC, y=log10p, color=sig)) +
  geom_point(data = plt_res[plt_res$sig == 'No', ], size=3) +
  geom_point(data = plt_res[plt_res$sig == 'Yes', ], size=3) +
  geom_point(data = plt_res[plt_res$sig == 'CredibleSet', ], size=3) +
  geom_text_repel(data = plt_res[plt_res$sig == 'CredibleSet', ], 
                  aes(label=variant_id), color='black',
                  force=1,  box.padding=unit(1,'lines'), 
  ) +
  scale_color_manual(values = c('CredibleSet'='#6a3d9a',
                                'Yes'='#cab2d6',
                                'No'='grey')) +
  theme_bw() +
  scale_x_continuous(limits = c(-2,2)) +
  labs(title=paste0('AD MPRA Locus: ',gwas_gene,' QTL: ',gene_name),
       x='logFC',
       y='-log10p'
  ) +
  theme(legend.position = 'none',
        axis.text = element_text(color='black', size=12),
        axis.title = element_text(color='black', size=12)
  )

plt_volcano

ggsave(filename, plot = plt_volcano, units = c("in"),
       width=8 ,height=8, dpi=600,
       limitsize = FALSE)










