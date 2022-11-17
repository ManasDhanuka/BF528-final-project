#------------------------------------------------------------------------------#

# Biological in-depth analysis

# Author: Manas Dhanuka

#------------------------------------------------------------------------------#

library(hgu133plus2.db)
library(tibble)
library(dplyr)
library(tidyverse)
library(GSEABase)

#importing file
de_matrix <- read.csv("project 1/n_biologist_data.csv", header = TRUE)

pset_ids <- de_matrix$probeids
pset_ids <- as.character(pset_ids)
mapped_results <- AnnotationDbi::select(hgu133plus2.db, keys=pset_ids, columns=("SYMBOL"))

new_de_matrix <- merge(x=de_matrix, y=mapped_results, by.x="probeids", by.y="PROBEID")
summary(new_de_matrix)

#removing rows with NA's
new_de_matrix <- new_de_matrix[!is.na(new_de_matrix$SYMBOL),]

gene_symbols <- new_de_matrix$SYMBOL
length(gene_symbols)
unique_gene_symbols <- unique(gene_symbols)
length(unique_gene_symbols)

# for gene symbol
new_de_matrix <- as_tibble(new_de_matrix)
removed_duplicates <- new_de_matrix %>%
  group_by(SYMBOL) %>%
  filter(p_adj == max(p_adj)) %>%
  ungroup(SYMBOL)

# Up-regulated
up <- slice_max(removed_duplicates, order_by=removed_duplicates$t_stat, n=1000)
up_regulated <- slice_max(removed_duplicates,order_by=removed_duplicates$t_stat, n=10)

# Down-regulated
down <- slice_min(removed_duplicates, order_by=removed_duplicates$t_stat, n=1000)
down_regulated <- slice_min(removed_duplicates,order_by=removed_duplicates$t_stat, n=10)

#writing to csv file
write.csv(up_regulated, "project 1/up_regulated.csv", row.names = FALSE)
write.csv(down_regulated, "project 1/down_regulated.csv", row.names = FALSE)

#get data from .gmt files
kegg <- getGmt("project 1/gmt_files/c2.cp.kegg.v7.2.symbols.gmt")
n_kegg <- length(kegg) # number of kegg gene sets
go <- getGmt("project 1/gmt_files/c5.all.v7.2.symbols.gmt")
n_go <- length(go) # number of go gene sets
hallmark <- getGmt("project 1/gmt_files/h.all.v7.2.symbols.gmt")
n_hallmark <- length(hallmark) # number of hallmark gene sets

#printing results
print("Number of gene sets in KEGG:")
print(n_kegg)
print("Number of gene sets in G0:")
print(n_go)
print("Number of gene sets in Hallmark:")
print(n_hallmark)

#filtering by p adjusted
DE <- removed_duplicates[removed_duplicates$p_adj < 0.05,]
DE_symbols <- DE$SYMBOL
all_gene_symbols <- removed_duplicates$SYMBOL

#symbols not in our DE data
non_DE_symbols <- setdiff(all_gene_symbols, DE_symbols)

#printing results
print("Total number of DE genes:")
length(DE_symbols)
print("Total number of non DE genes:")
length(non_DE_symbols)
print("Total number of genes left after removing duplicates:")
length(all_gene_symbols)

#function to create contingency table
contingency_function <- function(gene_set, de_gene_set, non_de_gene_set){
  
  # total number of differentially expressed genes
  n_differentially_expressed <- length(geneIds(de_gene_set))
  
  # total number of not differentially expressed genes
  n_not_differentially_expressed <- length(geneIds(non_de_gene_set))
  
  # number of differentially expressed genes in the gene set
  differentially_expressed_in <- gene_set & de_gene_set
  n_differentially_expressed_in <- length(geneIds(differentially_expressed_in))
  
  # number of differentially expressed genes not in the gene set
  n_differentially_expressed_not_in <- n_differentially_expressed - n_differentially_expressed_in
  
  # number of not differentially expressed genes in the gene set
  not_differentially_expressed_in <- gene_set & non_de_gene_set
  n_not_differentially_expressed_in <- length(geneIds(not_differentially_expressed_in))
  
  # number of not differentially expressed genes not in the gene set
  n_not_differentially_expressed_not_in <- n_not_differentially_expressed - n_not_differentially_expressed_in
  
  #concating data
  contingency_table <- c(n_differentially_expressed_in, n_differentially_expressed_not_in, n_not_differentially_expressed_in, n_not_differentially_expressed_not_in)
  
  return(contingency_table)
}

#creating gene set
gene_set_de <- GeneSet(DE_symbols, setName = "differentially expressed")
gene_set_not_de <- GeneSet(non_DE_symbols, setName = "not differentially expressed")

create_contingency <- function(gmt_data, gene_set_de, gene_set_not_de){
  
  contingency_gs <- list()
  for (i in gmt_data) {
    
    contigency <- contingency_function(i, gene_set_de, gene_set_not_de)
    
    contingency_gs[[setName(i)]] <- contigency
  }
  
  return(contingency_gs)
}

contingency_kegg <- create_contingency(kegg, gene_set_de, gene_set_not_de)
contingency_go <- create_contingency(go, gene_set_de, gene_set_not_de)
contingency_hallmark <- create_contingency(hallmark, gene_set_de, gene_set_not_de)

#for performing fisher test
fisher_ttest <-function(contingency_gs){
  
  gs_fisher_test <- data.frame(gene_set = character(), 
                               statistic_estimate = double(), 
                               p_value = double())
  
  for (i in names(contingency_gs)) {
    
    fisher_test <- fisher.test(matrix(contingency_gs[[i]], nrow = 2))
    gene_set = i
    statistic_estimate <- fisher_test$estimate[[1]]
    p_value <- fisher_test$p.value
    gs_fisher_test[nrow(gs_fisher_test) + 1,] <- c(i, statistic_estimate, p_value)
  }
  
  gene_set_name <- names(contingency_gs)
  gs_fisher_test$gene_set <- gene_set_name
  
  #adjusting p-value 
  gs_fisher_test$pvalue_adjusted <- p.adjust(gs_fisher_test$p_value, method = 'BH')
  
  #filtering significant values
  gs_significant <- gs_fisher_test[gs_fisher_test$p_value<0.05,]
  
  return(gs_significant)
  
}

#performing the fisher t_stat-test
kegg_fisher_test <- fisher_ttest(contingency_kegg)
go_fisher_test <- fisher_ttest(contingency_go)
hallmark_fisher_test <- fisher_ttest(contingency_hallmark)

print("Number of significantly enriched gene sets in KEGG:")
# print(n_kegg_significant)
print(length(kegg_fisher_test$gene_set))
print("Number of significantly enriched gene sets in GO:")
# print(n_go_significant)
print(length(go_fisher_test$gene_set))
print("Number of significantly enriched gene sets in Hallmark:")
# print(n_hallmark)
print(length(hallmark_fisher_test$gene_set))

print("Total number of significantly enriched genes:")
# print(n_kegg_significant + n_go_significant + n_hallmark_significant)
print(length(kegg_fisher_test$gene_set)+length(go_fisher_test$gene_set)+length(hallmark_fisher_test$gene_set))

#splicing top gene sets
top3_kegg <- slice_min(kegg_fisher_test, order_by = p_value, n=3)
top3_go <- slice_min(go_fisher_test, order_by = p_value, n=3)
top3_hallmark <- slice_min(hallmark_fisher_test, order_by = p_value, n=3)

top3_results <- rbind(top3_kegg, top3_go, top3_hallmark)
top3_results
write.csv(top3_results, "project 1/top3_results.csv", row.names = FALSE)

