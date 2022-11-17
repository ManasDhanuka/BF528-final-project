#------------------------------------------------------------------------------#

# Identifying differentially expressed genes associated with 
# myocyte differentiation 
# Author: Manas Dhanuka

#------------------------------------------------------------------------------#

library(dplyr)
library(tibble)
library(readr)
library(tidyverse)

file_path <- "/projectnb/bf528/users/im_not_dead_yet/project-2-im-not-dead-yet/P0_1_tophat/cuffdiff_out/gene_exp.diff"

#read in data
cuffdiff_data <- read.table(file_path, header = TRUE) 

#sort data by q-value
ordered_cuffdiff_data <- cuffdiff_data[order(cuffdiff_data[,"q_value"]),]

top10_genes <- ordered_cuffdiff_data[1: 10,]

#selecting names, FPKM values, log fold change, p-value, and q-value
top10_genes <- dplyr::select(top10_genes, gene, value_1, value_2, log2.fold_change., p_value, q_value)


hist(cuffdiff_data$log2.fold_change.,breaks = 30, xlab = "Log2 fold change",
     main = "Histogram of Log2_fold_change for all genes")
nrow(cuffdiff_data)


sig_data <- subset(cuffdiff_data, cuffdiff_data$significant == "yes")
nrow(sig_data)

P_Data <- subset(cuffdiff_data, cuffdiff_data$p_value < 0.01)
nrow(P_Data)

up_sub_p <- subset(P_Data,log2.fold_change.>0)
dw_sub_p <- subset(P_Data,log2.fold_change.<0)

nrow(up_sub_p)
nrow(dw_sub_p)

#write.csv(up_sub$gene,"up-regulated_gene.csv")
#write.csv(dw_sub$gene,"dw-regulated_gene.csv")

hist(sig_data$log2.fold_change.,breaks = 30, xlab = "Log2 fold change",
     main = "Histogram of Log2_fold_change for significant genes")

up_sub <- subset(sig_data,log2.fold_change.>0)
dw_sub <- subset(sig_data,log2.fold_change.<0)

#upregulated genes
nrow(up_sub)
#downregulated genes
nrow(dw_sub)

write.csv(up_sub$gene,"project 2/up-regulated_gene.csv")
write.csv(dw_sub$gene,"project 2/down-regulated_gene.csv")