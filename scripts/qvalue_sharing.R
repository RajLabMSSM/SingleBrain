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


### Plot for a specific result
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

# beta correlation and Rb values

#' Function for Rb analysis
#'
#' param b1 Beta from first dataset.
#' param se1 Standard error of beta from first dataset.
#' param b2 Beta from second dataset.
#' param se2 Standard error of beta from second dataset.
#' param theta Variable representing sample overlap between two datasets. Should be set 0 if no sample overlap.
#'
#' return Data frame with Rb, SE(Rb) and corresponding P-value.
#' export
#'
#' note This function comes from https://github.com/molgenis/sc-MetaBrain-consortium.
#'
#' examples
calcu_cor_true <- function(b1, se1, b2, se2, theta) {
  idx <- which(is.infinite(b1) | is.infinite(b2) | is.infinite(se1) | is.infinite(se2))
  if (length(idx) > 0) {
    b1 <- b1[-idx]
    se1 <- se1[-idx]
    b2 <- b2[-idx]
    se2 <- se2[-idx]
    theta <- theta[-idx]
  }

  var_b1 <- var(b1, na.rm = T) - mean(se1^2, na.rm = T)
  var_b2 <- var(b2, na.rm = T) - mean(se2^2, na.rm = T)
  if (var_b1 < 0) {
    var_b1 <- var(b1, na.rm = T)
  }
  if (var_b2 < 0) {
    var_b2 <- var(b2, na.rm = T)
  }
  cov_b1_b2 <- cov(b1, b2, use = "complete.obs") - mean(theta, na.rm = T) * sqrt(mean(se1^2, na.rm = T) * mean(se2^2, na.rm = T))
  r <- cov_b1_b2 / sqrt(var_b1 * var_b2)

  r_jack <- c()
  n <- length(b1)
  for (k in 1:n) {
    b1_jack <- b1[-k]
    se1_jack <- se1[-k]
    var_b1_jack <- var(b1_jack, na.rm = T) - mean(se1_jack^2, na.rm = T)
    b2_jack <- b2[-k]
    se2_jack <- se2[-k]
    var_b2_jack <- var(b2_jack, na.rm = T) - mean(se2_jack^2, na.rm = T)
    if (var_b1_jack < 0) {
      var_b1_jack <- var(b1_jack, na.rm = T)
    }
    if (var_b2_jack < 0) {
      var_b2_jack <- var(b2_jack, na.rm = T)
    }
    theta_jack <- theta[-k]
    cov_e1_jack_e2_jack <- mean(theta_jack, na.rm = T) * sqrt(mean(se1_jack^2, na.rm = T) * mean(se2_jack^2, na.rm = T))
    cov_b1_b2_jack <- cov(b1_jack, b2_jack, use = "complete.obs") - cov_e1_jack_e2_jack
    r_tmp <- cov_b1_b2_jack / sqrt(var_b1_jack * var_b2_jack)
    r_jack <- c(r_jack, r_tmp)
  }
  r_mean <- mean(r_jack, na.rm = T)
  idx <- which(is.na(r_jack))
  if (length(idx) > 0) {
    se_r <- sqrt((n - 1) / n * sum((r_jack[-idx] - r_mean)^2))
  } else {
    se_r <- sqrt((n - 1) / n * sum((r_jack - r_mean)^2))
  }

  p <- pchisq((r / se_r)**2, df = 1, lower.tail = FALSE)

  res <- cbind(r, se_r, p)
  return(res)
}

rb_res_all <- data.frame()
rb_se_p_all <- data.frame()

### Set colorset on brain celltype
cell_list <- c('Ast','End',
               'Ext','IN','MG',
               'OD','OPC'
                      )

colorset <- c('#1F78B4',    # Ast
            '#b15928',    # End
            '#33A02C',    # Ext
            '#E31A1C',    # IN
            '#6a3d9a',    #MG
            '#FF7F00',    # OD
            '#FDBF6F'    # OPC
            )
names(colorset) = cell_list

## Collate data from qvalue_sharing RData (Take a look: https://github.com/RajLabMSSM/downstream-QTL).
for(celltype in cell_list){
  message(celltype)
  load(paste0(path_dir,"SingleBrain_",celltype,"_expression_EUR:",
              "Fujita_",celltype,"_eQTL:0.05.qvalue.RData"))

  target_res <- target_res[(target_res$target_beta!=0)&(target_res$target_pvalue!=1),]

  rb_res <- merge(source_top_res, target_res)
  rb_res <- rb_res %>% dplyr::filter((source_pvalue!=0)&(target_pvalue!=0)) %>% 
    mutate(sig = case_when(
    target_p_perm < 0.05 ~ 'sig',
    target_p_perm >==  0.05 ~ 'no'
    )) %>% mutate(target_beta_fixed = case_when(
    (source_a1==target_a1)&(source_a2==target_a2)  ~ target_beta,
    (source_a1!=target_a1)&(source_a2!=target_a2) ~ -target_beta
    ))

  source_max <- (abs(rb_res$source_beta))
  target_max <- (abs(rb_res$target_beta))

  rb_se_p <- calcu_cor_true(rb_res$source_beta, rb_res$source_se,
                            rb_res$target_beta, rb_res$target_se, 
                          overlap_samples) %>% as.data.frame()
  
  rb_se_p$nsnp <- nrow(rb_res)
  rb_se_p$celltype <- celltype
  
  rb_res$celltype <- celltype

  message(res_cor$estimate)
  rb_res <- rb_res %>% dplyr::filter(sig=='sig')
  message(nrow(rb_res))
  rb_res <- sig_filter_singlebrain(celltype, rb_res)
  message(nrow(rb_res))
  
  rb_res_all <- rbind(rb_res_all, rb_res)
  rb_se_p_all <- rbind(rb_se_p_all, rb_se_p)
  
}
rb_res_all$celltype <- factor(rb_res_all$celltype, levels = c('Ext','IN','OD','OPC',
'Ast','MG','End'))

rb_se_p_all # check Rb values 

### Beta correlation plot between two cohorts by each celltype
x_y_max <- max(abs(rb_res_all$source_beta), abs(rb_res_all$target_beta))

p_beta <- rb_res_all %>% ggplot(aes(x=source_beta, y=target_beta_fixed, color=celltype)) + 
  geom_point(size=1, alpha=0.5) +
  scale_color_manual(values = colorset) +
  scale_x_continuous(limits = c(-x_y_max, x_y_max)) +
  scale_y_continuous(limits = c(-x_y_max, x_y_max)) +
  facet_wrap(. ~ celltype, ncol = 4
                   ) +
  labs(x = "beta of SingleBrain", y = "beta of Fujita et al.") +
  theme_bw() +
  theme(title = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, colour = "black"), 
        axis.title.x = element_text(size=12, colour = "black"),
        axis.text = element_text(size=12, colour = "black"), 
        legend.position = 'none' ,
        panel.spacing = unit(5, "mm"),
        strip.placement = "outside",
        strip.text = element_text(colour = "black",  size=12),
        strip.background = element_blank(),
        )
  
show(p_beta)
filename <- paste0('~/beta_singlebrain_fujita.pdf')
gsave(filename, plot = p_beta)

