suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(edgeR))

option_list <- list(
  make_option(c('-c', '--cell'), help='the celltype', default = "")
)

  option.parser <- OptionParser(option_list=option_list)
  opt <- parse_args(option.parser)
  
  celltype <- opt$cell
  print(celltype)

  print("Psuedobulk processing")

  dataset <- 'SingleBrain'
  path_dir='~/'
  assay_name <- 'RNA'

  df_meta <- fread(file=paste0(path_dir, celltype, "/metadata.csv"))  %>% as.data.frame()
  row.names(df_meta) <- df_meta$V1
 
  getwd()
  data_10x <- ReadMtx(
    mtx =  paste0(path_dir, celltype, '/matrix.mtx.gz'),
    features = paste0(path_dir, celltype, "/features.tsv.gz"),
    cells = paste0(path_dir, celltype, "/barcodes.tsv.gz")
  )
  data.seurat <- CreateSeuratObject(counts = data_10x, meta.data = df_meta)

  ### Filter individual under 10 cells
  tbl <- table(data.seurat@meta.data$sampleID)
  df_counts <- as.data.frame(tbl)
  df_filter <- df_counts %>% filter( Freq >= 10 ) %>% select(Var1)
  colnames(df_filter) <- c('sample_id')
   
  sample_list <- unique(data.seurat@meta.data[!(is.na(data.seurat@meta.data$sampleID)),'sampleID' ])
  gene_key <- '~/gene_id_key.txt'
  df_gene <- fread(gene_key, sep=' ')
  colnames(df_gene) <- c('gene_id', 'transcript_id')
  
  print("Making Mean")
  dim(data.seurat)
  
  mean.cell <-
    purrr::map( sample_list, ~{
      tmp <- subset(data.seurat, sampleID == .x)
      mean_tmp <- data.frame(rowMeans(GetAssayData(tmp, assay = assay_name ,slot = "counts")))
      colnames(mean_tmp) <- .x
      return(mean_tmp)
    })
  mean.cell <- as.data.frame(mean.cell)
  mean.cell$transcript_id <- rownames(mean.cell)
  print(nrow(mean.cell))
  df_result <- merge(df_gene, mean.cell, by='transcript_id' )[,-1]
  
  name = df_result[,1]$gene_id
  
  genes_tpm <- data.frame(df_result, row.names = name )
  
  row.names(genes_tpm) <- NULL
  print(nrow(genes_tpm))
  
  path_result <- paste0(path_dir, dataset, "_" ,celltype,"_",assay_name,"_mean.tsv")
  print(path_result)
  write.table(genes_tpm, path_result, row.names = F, quote = F, sep='\t')
  
  print("Making Sum")
  sum.cell <-
    purrr::map( sample_list, ~{
      tmp <- subset(data.seurat, sampleID == .x)
      sum_tmp <- data.frame(rowSums(GetAssayData(tmp, assay = assay_name, slot = "counts")))
      colnames(sum_tmp) <- .x
      return(sum_tmp)
    })
  sum.cell <- as.data.frame(sum.cell)
  sum.cell$transcript_id <- rownames(sum.cell)
  print(nrow(sum.cell))
  
  df_result <- merge(df_gene, sum.cell, by='transcript_id' )[,-1]
  name = df_result[,1]$gene_id
  
  genes_counts <- data.frame(df_result, row.names = name )
  
  row.names(genes_counts) <- NULL
  
  print(nrow(genes_counts))
  
  path_result <- paste0(path_dir, dataset, "_" ,celltype,"_",assay_name,"_counts.tsv")
  print(path_result)
  write.table(genes_counts, path_result, row.names = F, quote = F, sep='\t')
  
  print("Making RData")

  df_key <- fread('~/total_sample_key.txt') %>%
    dplyr::select(sample_id)
  df_key <- merge(df_key, df_filter)

  ### Sample key 
  df_sample_key <- fread('~/total_sample_key.txt')
  df_sample_key <- merge(df_sample_key, df_key, by <- 'sample_id')

  result_key <- paste0('~/',dataset,'_sample_key_', celltype,'.txt')
  print(result_key)
  write.table(df_sample_key, result_key, quote = F, row.names = F, sep='\t')
  

