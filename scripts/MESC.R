suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(readxl))

all_mesc <- read_tsv(paste0('~/singleBrain_ALL_MESC.tsv'))
all_mesc$cell_type <- factor(levels=c('Ext','IN','Ast','OD','OPC','MG','End'), all_mesc$cell_type)

### Set color set
dataset_celltype <- c('Ast','MG',
                      'Ext','IN',
                      'OD','OPC',
                      'End')

colorset <- c('#1F78B4',    # Astrocyte
              '#6a3d9a',    # Microglia
              '#33A02C',    # Excitatory
              '#E31A1C',    # Inhibitory
              '#FF7F00',    # Oligodendrocyte
              '#FDBF6F',    # OPC
              '#b15928'    # Endothelial
)
names(colorset) = dataset_celltype

colorset

## MESC plot 
p_mesc <- all_mesc %>% dplyr::filter( disease %in% c("SCZ", "MS", "AD",'PD','BPD', 'ALS')) %>% dplyr::mutate(disease=factor(disease, levels=c("SCZ", "MS", "AD",'PD','BPD', 'ALS'))) %>% 
  ggplot(aes(x=cell_type)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2,) +
  geom_errorbar(aes(ymin = h2med_h2g - SE, ymax = h2med_h2g + SE), width = 0.25, size = 0.5 ) + 
  geom_point(aes(x = cell_type, y = h2med_h2g, colour = cell_type), size = 2  ) + 
  scale_colour_manual(values = colorset) +
  facet_grid( ~ disease, scales = "free", space = "free_y") + 
  theme_classic() + 
  coord_flip() +
  theme(axis.title=element_text(size=10, color='black'),
        axis.text.x=element_text(size=10, color='black',angle = 90),
        axis.text.y=element_text(size=10, color='black'),
        axis.ticks.x = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        legend.ticks = element_blank(),
        legend.position = 'none',
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(),
        strip.text = element_text(color = "black", size=12),
        axis.text = element_text(colour = "black"),
        plot.background = element_rect(fill="white", colour='NA'), 
  ) +
  labs(y = expression(italic(h)["med"]^2 / italic(h)["g"]^2 ), x = "" ) +
  scale_x_discrete(limits = rev) +
  theme(strip.text.y = element_text(angle = 0))

show(p_mesc)


filename <- paste0('~/errorbarplot_mesc_eqtl_main.pdf')
print(filename)
ggsave(filename, plot = p_mesc, units = c("in"),
       width=7 ,height=3, dpi=600,
       limitsize = FALSE)
