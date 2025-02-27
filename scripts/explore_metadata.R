library(tidyverse)
library(ggplot2)
library(splitstackshape)
library(data.table)
library(dplyr)


## Circle figure

data <-read_tsv('~/cell_count.csv')

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nObsType <- nlevels(as.factor(data$celltype))
# to_add <- data.frame( matrix(NA, empty_bar*nlevels(df_counts$cohort)*nObsType, ncol(data)) )
to_add <- data.frame( matrix(NA, empty_bar*4*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$cohort <- rep(levels(data$cohort), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(cohort, individualID)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data <- data %>% group_by(id, individualID) %>% summarize(tot=sum(Freq))
number_of_bar <- nrow(label_data)
# I substract 0.5 because the letter must have the angle of the center of the bars. 
# Not extreme right(1) or extreme left (0)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(cohort) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end))) 

base_data2 <- data %>% 
  group_by(cohort) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end))) %>% na.omit()


# prepare a data frame for grid (scales)
grid_data <- base_data2
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

p_circle <- ggplot(data) +      
  # Add the stacked bar
  # geom_bar(aes(x=as.factor(id), y=Freq, fill=celltype), stat="identity") +
  geom_bar(aes(x=as.factor(id), y=Freq, fill=celltype), stat="identity") +
  scale_fill_manual(values = colorset) + 
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), 
               colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), 
               colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), 
               colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 150, xend = start, yend = 150), 
               colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), 
               = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  
  ylim(-2, max(label_data$tot, na.rm=T)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add base line information
  geom_segment(data=base_data2, aes(x = start, y = -0.1, xend = end, yend = -0.1),
               colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data2, aes(x = title, y = -0.2, label=cohort), hjust=c(1,0,0,0), 
            colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  theme(legend.position = 'right')

p_circle

filename <- paste0(path_dir, '/figures/','singlebrain_circle_stacked_bar.pdf')
print(filename)
ggsave(filename, plot = p_circle, width=30 ,height=30, dpi=600)


# Metadata summary
df_meta <- read_tsv('~/metadata.csv')

## Sex
dataset_list <- c('Fujita','Mathys','Gabitto','Bryois')

df_tmp$dataset <- factor(df_tmp$dataset, levels = dataset_list)
# df_tmp 
df_tmp <- df_tmp %>% group_by(dataset, sex) %>%
  count(sex) 
df_tmp <- df_tmp  %>% group_by(dataset) %>% mutate(frac = n / sum(n) )
df_tmp$pct <- format(round(df_tmp$frac*100, 1), nsmall = 1) 
df_tmp$labels <- paste0(df_tmp$n, '\n(',df_tmp$pct,'%)') 

p_sex <- ggplot(df_tmp, aes(dataset, frac*100, fill = sex)) +
  geom_col() +
  theme_classic() +
  scale_fill_manual(breaks = c('Male','Female'), values = c('#1f78b4','#e31a1c')) +
  labs(x='Dataset', y='Percentage of sex (%)') +
  theme(axis.text = element_text(color='black', size=10),
        axis.text.x = element_text(face = 'italic'),
        axis.title = element_text(color='black', size=12),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(color='black', size=10),
  ) +
  geom_text(
    aes(label = labels),
    position=position_stack( vjust = 0.5),
    color='black' , size=3
  ) 

filename <- paste0('~/supple_sex_stackbar.pdf')
ggsave(filename, plot = p_sex, width=4.5 ,height=5, dpi=600)

## Age
df_tmp <- df_meta %>% dplyr::select(dataset, sample_id, ageDeath, sex)
p_age <- ggplot(df_tmp, aes(x=dataset, y=ageDeath, fill=sex)) + 
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(breaks = c('Male','Female'), 
                    values = c('#1f78b4','#e31a1c')) +
  labs(x='Dataset', y='Age at death') +
  theme(axis.text = element_text(color='black', size=10),
        axis.text.x = element_text(face = 'italic'),
        axis.title = element_text(color='black', size=12),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(color='black', size=10),
  ) 
p_age
filename <- paste0('~/supple_age_boxplot.pdf')
ggsave(filename, plot = p_age, width=4.5 ,height=5, dpi=600)

# Age2
p_agehisto <- ggplot(df_tmp, aes(x=ageDeath, fill=sex)) + 
  geom_histogram() +
  theme_bw() +
  scale_fill_manual(breaks = c('Male','Female'), 
                    values = c('#1f78b4','#e31a1c')) +
  labs(x='Age at death', y='Counts') +
  theme(axis.text = element_text(color='black', size=10),
        axis.title = element_text(color='black', size=12),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(color='black', size=10),
  ) 
p_agehisto
filename <- paste0('~/supple_age_histogram.pdf')
ggsave(filename, plot = p_agehisto, width=5 ,height=5, dpi=600)

## Case n Control
df_tmp <- df_meta %>% dplyr::select(dataset, sample_id, case_or_control)

p_case <- ggplot(df_tmp, aes(dataset, frac*100, fill = case_or_control)) +
  geom_col() +
  theme_classic() +
  scale_fill_manual(breaks = c('Control','Case'), 
                    values = c('#33a02c','#ff7f00')) +
  labs(x='Dataset', y='Percentage of sex (%)') +
  theme(axis.text = element_text(color='black', size=10),
        axis.title = element_text(color='black', size=12),
        axis.text.x = element_text(face = 'italic'),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(color='black', size=10),
  ) +
  geom_text(
    aes(label = labels),
    position=position_stack( vjust = 0.5),
    color='black' , size=3
  ) 

show(p_case)
filename <- paste0('~/supple_case_stackbar.pdf')
ggsave(filename, plot = p_case, width=4.5 ,height=5, dpi=600)

## Diagnosis
df_tmp <- df_meta %>% dplyr::select(dataset, sample_id, Diagnosis)

df_tmp <- df_tmp %>% group_by(Diagnosis) %>%
  count(Diagnosis) 
df_tmp$frac <- df_tmp$n / sum(df_tmp$n)
df_tmp$pct <- format(round(df_tmp$frac*100, 1), nsmall = 1) 
df_tmp$labels <- paste0(df_tmp$n, '\n(',df_tmp$pct,'%)') 

df_tmp$Diagnosis <- factor(df_tmp$Diagnosis, levels = c('Control','AD','MCI','Dementia','Other/Unknown'))

p_diagnosis <- ggplot(df_tmp, aes(Diagnosis, n, fill = Diagnosis)) +
  geom_col() +
  theme_classic() +
  scale_fill_manual(breaks = c('Control','AD','MCI','Dementia','Other/Unknown'), 
                    values = c('#33a02c','#e31a1c','#ff7f00','#1f78b4','#6a3d9a')
  ) +
  labs(x='Diagnosis', y='Diagnosis counts (n)') +
  theme(axis.text = element_text(color='black', size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.title = element_text(color='black', size=12),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(color='black', size=10),
  ) +
  geom_text(
    aes(label = labels),
    position=position_stack( vjust = 0.5),
    color='black' , size=3
  )

show(p_diagnosis)
filename <- paste0('~/supple_diagnosis_bar.pdf')
ggsave(filename, plot = p_diagnosis, width=4.5 ,height=5, dpi=600)

## Diagnosis stack bar
df_tmp <- df_meta %>% dplyr::select(dataset, sample_id, Diagnosis)
df_tmp <- df_tmp %>% group_by(dataset, Diagnosis) %>%
  count(Diagnosis) 
df_tmp <- df_tmp  %>% group_by(dataset) %>% mutate(frac = n / sum(n) )
df_tmp$pct <- format(round(df_tmp$frac*100, 1), nsmall = 1) 
df_tmp$labels <- paste0(df_tmp$n, ' (',df_tmp$pct,'%)') 

df_tmp$Diagnosis <- factor(df_tmp$Diagnosis, levels = c('Control','AD','MCI','Dementia','Other/Unknown'))
df_tmp$dataset <- factor(df_tmp$dataset, levels = dataset_list)

p_diagnosis <- ggplot(df_tmp, aes(dataset, frac*100, fill = Diagnosis)) +
  geom_col() +
  theme_classic() +
  scale_fill_manual(breaks = c('Control','AD','MCI','Dementia','Other/Unknown'), 
                    values = c('#33a02c','#e31a1c','#ff7f00','#1f78b4','#6a3d9a')
  ) +
  labs(x='Dataset', y='Percentage of diagnosis (%)') +
  theme(axis.text = element_text(color='black', size=10),
        axis.title = element_text(color='black', size=12),
        axis.text.x = element_text(face = 'italic'),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(color='black', size=10),
  ) +
  geom_text(
    aes(label = labels),
    position=position_stack( vjust = 0.5),
    color='black' , size=3
  )

filename <- paste0('~/supple_diagnosis_stackbar.pdf')
ggsave(filename, plot = p_diagnosis, width=5 ,height=7, dpi=600)


# PCA plot

df_pca <- read_tsv('~/somalier-ancestry.somalier-ancestry.tsv')
df_pca$participant_id <- df_pca$`#sample_id`

df_pca <- select(df_pca, participant_id, predicted_ancestry , EUR_prob, starts_with("PC"))

df_pca_1kg <- df_pca[1:2504,] 
df_pca_1kg$cohort <- '1kG'
df_pca_1kg$alphas <- 0.01
df_pca_1kg$sizes <- 2
df_pca_1kg$sample_id <- df_pca_1kg$participant_id
df_pca_1kg$sex <- NA
df_pca_1kg$case_or_control <- NA
df_pca_1kg$dataset <- '1kG'

df_pca_singlebrain <- df_pca[2505:nrow(df_pca),] 
df_pca_singlebrain$cohort <- 'SingleBrain'
df_pca_singlebrain$alphas <- 0.6
df_pca_singlebrain$sizes <- 3

all_res_pca <- df_pca_singlebrain

all_res_pca$Ancestry <- all_res_pca$predicted_ancestry

plot_pca <- all_res_pca %>% ggplot(aes(x=PC1, y=PC2, color=Ancestry, shape=cohort)) +
  geom_point(aes(size=sizes, alpha=alphas)) +
  scale_shape_manual(values = c(16, 15))+
  scale_color_manual(values = colorset) + 
  scale_alpha_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  scale_size_continuous(limits = c(1,5), breaks = c(1,2,3,4,5)) +
  theme_bw() +
  theme(axis.text.x = element_text(color='black',  size = 14),
        axis.text.y = element_text(color='black',  size = 14),
        axis.title.x = element_text(color='black',  size = 14),
        axis.title.y = element_text(color='black',  size = 14)
  )

show(plot_pca)

filename <- paste0('~/dotplot_pca_ancestry.pdf')
ggsave(filename, plot = plot_pca, width=7 ,height=6, dpi=600)











