#------------------------------------------------------------------------------#

# Data Pre-processing and Quality Control

# Author: Manas Dhanuka

#------------------------------------------------------------------------------#

library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)
library(ggfortify)
library(dplyr)

filepath = '/projectnb/bf528/users/im_not_dead_yet/project_1/samples'

celpath <- system.file('celfiles', package='affydata')
fns <- list.celfiles(path=celpath, full.names=TRUE)

# read in files
celbatch <- ReadAffy(celfile.path=filepath)

# normalize files together
RMA <- rma(celbatch)

# convert AffyBatch to PLMset
pset <- fitPLM(celbatch, normalize = TRUE, background = TRUE)

# relative log expression (RLE)
rle_stats <- data.frame(t(affyPLM::RLE(pset, type='stats')))

# plot rle_stats
png("RLE_hist.png")
rle_medians <- ggplot(rle_stats, aes(x=median)) + 
  geom_histogram(bins=50, 
                 color = 'black', 
                 fill = 'lightblue') +
  labs(title = 'Histogram of RLE Median Distribution') +
  theme_linedraw()
dev.off()

# normalized unscaled standard errors (NUSE)
nuse_stats <- data.frame(t(NUSE(pset, type = 'stats')))

# plot nuse_stats
png("NUSE_hist.png")
nuse_medians <- ggplot(nuse_stats, aes(x=median)) + 
  geom_histogram(bins=50, 
                 color = 'black', 
                 fill = 'lightblue') +
  labs(title = 'Histogram of NUSE Median Distribution') +
  theme_linedraw()
dev.off()

# correction for batch effects

#importing provided metadata
annotation_data <- read.csv(file = '/project/bf528/project_1/doc/proj_metadata.csv')

edata <- exprs(RMA)
batch <- annotation_data$normalizationcombatbatch

modcombat <- model.matrix(~as.factor(normalizationcombatmod), data = annotation_data)

combat_edata = ComBat(dat = edata, batch = batch, mod = modcombat)

# write combat_edata to csv
write.csv(combat_edata, file = 'project 1/input_data/edata.csv')

# transpose before pca
trans_edata <- t(combat_edata)

# scale and center
scaled_edata <- scale(trans_edata, center = TRUE, scale = TRUE)

# retranspose
scaled_edata <- t(scaled_edata)

# perform pca
pca <- prcomp(scaled_edata, center = FALSE, scale = FALSE)

# percent variance explained
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
var_explained[1:5]

# plot pca

png("pca_plot.png")

pca_plot <- pca$rotation %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(size = 1, color='blue') +
  theme_linedraw() +
  labs(x = paste0('PC1: ', round(var_explained[1]*100, 2), '%'),
       y = paste0('PC2: ', round(var_explained[2]*100, 2), '%'),
       title = 'PC1 vs PC2')

dev.off()

pca_plot