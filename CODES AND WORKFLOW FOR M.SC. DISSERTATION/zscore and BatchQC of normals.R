#Title: zscore and BatchQC of normals 
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

samples <- list.files(path = "C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/", full.names = T, pattern = "_quant")
files <- file.path(samples, "quant.sf")

names(files) <- str_replace(samples, "C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/", "") %>%
  str_replace("_quant","")

tx2gene <- read.delim("C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/txt2gene_ENST_genesymbol.csv", sep = ",", col.names = c("Transcript_id", "Gene_symbol", "Gene_id"))
tx2gene %>% View()

txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("Transcript_id", "Gene_symbol", "Gene_id")], countsFromAbundance="lengthScaledTPM")

# Count matrix generation
counts <- txi$counts %>% 
  round() %>% 
  data.frame()

#countsdata file
write.csv(counts, file = "C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/HNSCC_counts_all_normals.csv", sep = "\t", col.names = TRUE, row.names = TRUE)

# Length matrix generation
length <- txi$length %>% 
  round() %>% 
  data.frame()


#lengthdata file
write.csv(length, file = "C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/HNSCC_length_all_normals.csv", sep = "\t", col.names = TRUE, row.names = TRUE)

#metadata
meta <-read.csv("C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/Meta_HNSCC_NORMALS.csv",row.names=1)
head(meta)
dim(meta)
colnames(meta)
rownames(meta)

all(rownames(meta)==colnames(counts))
meta$Status <- factor(meta$Status)

data <- read.csv("C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/HNSCC_zscore_median.csv")
median <- data$Median

combat = ComBat_seq(as.matrix(counts), batch = meta$Status)
write.csv(combat, file = "C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/HNSCC_combat_20230824.csv", sep = ",")

combat <- as.matrix(combat)

dim(combat)

tpm_mat <-tpm(combat,median)
write.csv(tpm_mat, file = "C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/HNSCC_tpm_matrix_20230824.csv", sep = "\t", col.names = TRUE, row.names = TRUE)

#z_score
Z <- t(scale(t(tpm_mat)))
Z_score <- write.csv(Z, file = "C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/z_score_matrix_20230824.csv", sep = "\t", col.names = TRUE, row.names = TRUE)

Batchqc = batchQC(counts, meta$Status,  
                  report_file="C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/batchqc_before_combat_report_20230827.html", report_dir=".",
                  report_option_binary="111111111",
                  view_report=FALSE, interactive=TRUE, batchqc_output=TRUE)

Batchqc_combat = batchQC(tpm_mat, meta$Status,  
                         report_file="C:/Users/amrit/Desktop/HNSCC_NORMALS_ANALYSIS_20230816/batchqc_after_combat_report_20230827.html", report_dir=".",
                         report_option_binary="111111111",
                         view_report=FALSE, interactive=TRUE, batchqc_output=TRUE)

boxplot(log2(counts + 1),xlab="Samples",ylab = "Counts",names.arg = colnames(data), beside = TRUE, cex.axis = 0.4, cex.label =0.6, las=2,ylim= c(-2,10))

boxplot(log2(tpm_mat + 1),xlab="Samples",ylab = "TPM",names.arg = colnames(data), beside = TRUE, cex.axis = 0.4, cex.label =0.6, las=2,ylim= c(-2,10))



