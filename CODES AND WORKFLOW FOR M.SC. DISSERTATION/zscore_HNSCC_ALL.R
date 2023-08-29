#Title: zscore of normals and tumors
#Date: 17-08-2023

remotes::install_github("davidrequena/drfun")
library(DRnaSeq)
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(stringr)
library(affy)
library(scales)
library(sva)
library(BatchQC)
library(graphics)

#The variables store the quant file corresponding to each sample, the same can be viewed in files

samples <- list.files(path = "C:/Users/amrit/Desktop/HNSCC_ALL/", full.names = T, pattern = "_quant")
files <- file.path(samples, "quant.sf")

names(files) <- str_replace(samples, "C:/Users/amrit/Desktop/HNSCC_ALL/", "") %>%
  str_replace("_quant","")

tx2gene <- read.delim("C:/Users/amrit/Desktop/HNSCC_ALL/txt2gene_ENST_genesymbol.csv", sep = ",", col.names = c("Transcript_id", "Gene_symbol", "Gene_id"))
tx2gene %>% View()

txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("Transcript_id", "Gene_symbol", "Gene_id")], countsFromAbundance="lengthScaledTPM")

# Count matrix generation
counts <- txi$counts %>% 
  round() %>% 
  data.frame()

#countsdata file
write.csv(counts, file = "C:/Users/amrit/Desktop/HNSCC_ALL/HNSCC_ALL_20230824.csv", sep = "\t", col.names = TRUE, row.names = TRUE)

length <- txi$length %>% 
  round() %>% 
  data.frame()

#lengthdata file
write.csv(length, file = "C:/Users/amrit/Desktop/HNSCC_ALL/HNSCC_length_20230824.csv", sep = "\t", col.names = TRUE, row.names = TRUE)

#metadata
meta <-read.csv("C:/Users/amrit/Desktop/HNSCC_ALL/Meta_HNSCC_ALL.csv",row.names=1)
head(meta)
dim(meta)
colnames(meta)
rownames(meta)

all(rownames(meta)==colnames(counts))
meta$Status <- factor(meta$Status)

data <- read.csv("C:/Users/amrit/Desktop/HNSCC_ALL/HNSCC_zscore_median.csv")
median <- data$Median

combat = ComBat_seq(as.matrix(counts), batch = meta$Status)
write.csv(combat, file = "C:/Users/amrit/Desktop/HNSCC_ALL/HNSCC_combat_20230824.csv", sep = ",")

combat <- as.matrix(combat)

dim(combat)

tpm_mat <-tpm(combat,median)
write.csv(tpm_mat, file = "C:/Users/amrit/Desktop/HNSCC_ALL/HNSCC_tpm_matrix_20230824.csv", sep = "\t", col.names = TRUE, row.names = TRUE)

#z_score
Z <- t(scale(t(tpm_mat)))
Z_score <- write.csv(Z, file = "C:/Users/amrit/Desktop/HNSCC_ALL/z_score_matrix_20230824.csv", sep = "\t", col.names = TRUE, row.names = TRUE)

