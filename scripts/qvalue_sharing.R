suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))


# qvaule sharing - within SingleBrain
major <- c('Ext', 'IN','Ast', 'OD',  'OPC', 'MG', 'End')
miga <- c('Ext', 'IN', 'Ast', 'OD',  'OPC', 'MG', 'End','MiGA')
qtl_list <- c('Ext', 'Ext1', 'Ext2', 'Ext4',  'Ext8', 'Ext7', 'Ext5', 'Ext3', 'Ext6',
              'IN',  'IN4', 'IN7', 'IN6', 'IN2','IN3','IN5','IN1',  
              
              'OD',  'OD4', 'OD2',  'OD1', 
              'Ast','Ast1', 'Ast2','Ast3', 'Ast4',
              'OPC',  'OPC1', 'OPC2', 
              'MG', 'MG1', 'MG2',  'MG3',  'MG4',
              'End')


qvalues<- read_tsv("~/qvalue_res_singlebrain/all_qvalue_merged.0.05.tsv")

qvalues_to_plot <-
  qvalues %>%
  mutate(pi1 = ifelse(is.na(pi1), yes = 1, no = pi1)) %>%
  mutate(pi1 = round(pi1, digits = 2)) %>%
  mutate(pi1 = ifelse( source_name == target_name, yes = 1, no = pi1))

qvalues_to_plot$pi1 <- round(qvalues_to_plot$pi1, digits = 2)

options(digits=2)
qvalues_to_plot_major <- qvalues_to_plot[which(qvalues_to_plot$source_name %in% major ),]
qvalues_to_plot_major <- qvalues_to_plot_major[which(qvalues_to_plot_major$target_name %in% major ),]

qvalues_to_plot_major$source_name <- factor(qvalues_to_plot_major$source_name , levels = major)
qvalues_to_plot_major$target_name <- factor(qvalues_to_plot_major$target_name , levels = major)
qvalues_to_plot_major$pi1 <- round(qvalues_to_plot_major$pi1, digits = 2)

plot_q <- qvalues_to_plot_major %>% 
  ggplot(aes(x= source_name, y = target_name, fill = pi1, 
             label=sprintf("%0.2f", round(pi1, digits = 2)))) + 
  geom_tile() + geom_text() +
  labs(x = "Discovery", y = "Replication") +
  scale_fill_gradientn(colours = c("#1f78b4","#ffff99","#e31a1c"), 
                       limits = c(0,1),
                       breaks = c(0,1),
                       labels = c(0,1),
  ) +
  theme_void() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # coord_flip() +
  scale_y_discrete(limits = rev(levels(qvalues_to_plot_major$target_name))) +
  theme(axis.title.y = element_text(size=12, colour = "black", angle=90), 
        axis.title.x = element_text(size=12, colour = "black", angle=0),
        axis.text = element_text(size=12, colour = "black"), 
        text = element_text(size=12, colour = "black"),
        axis.ticks = element_blank()
  )
show(plot_q)

ggsave(filename, plot = plot_q, 
       width=5 ,height=4, 
       dpi=600)


qvalues_to_plot$source_name <- factor(qvalues_to_plot$source_name , levels = qtl_list)
qvalues_to_plot$target_name <- factor(qvalues_to_plot$target_name , levels = qtl_list)


plot_q <- qvalues_to_plot %>%
  ggplot(aes(x= source_name, y = target_name, 
             fill = pi1, 
             label=sprintf("%0.2f", round(pi1, digits = 2)))) + 
  geom_tile() + 
  geom_text() +
  labs(x = "Discovery", y = "Replication") +
  scale_fill_gradientn(colours = c("#1f78b4","#ffff99","#e31a1c"),
                       # scale_fill_distiller(palette = "Spectral",
                       limits = c(0,1),
                       breaks = c(0,1),
                       labels = c(0,1),) +
  theme_void() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # coord_flip() +
  scale_y_discrete(limits = rev(levels(qvalues_to_plot$target_name))) +
  # scale_x_discrete(expand = c(0,0)) + scale_y_discrete( expand = c(0,0)) + 
  theme(axis.title.y = element_text(size=12, colour = "black", angle=90), 
        axis.title.x = element_text(size=12, colour = "black", angle=0),
        axis.text = element_text(size=12, colour = "black"), 
        text = element_text(size=8, colour = "black"),
        axis.ticks = element_blank()
  )

show(plot_q)

ggsave(filename, plot = plot_q, 
       width=16 ,height=15, 
       dpi=600)





