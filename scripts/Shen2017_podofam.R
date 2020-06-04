#!/usr/bin/env Rscript

### Measuring phylogenetic signal of markers for the Podosporaceae family
#############################################################################

# The aim of this script is to plot the values of Site and Gene-wise
# log-likelihood scores (SLS and GLS) of a number of markers

# Part of the Podosporaceae.smk pipeline

# Based on the scripts of Iker Irisarri, which in turn takes output from the
# scripts in Shen et al. (2017) Contentious relationships in phylogenomic
# studies can be driven by a handful of genes, Nat. Ecol. Evol

# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020-06-04
# Version 1
# =======================================

library(ggplot2, quietly = TRUE)
library(cowplot) 
library(tidyr) # for gather

# --- Data
podo_SLSdata_gene <- snakemake@input$sitwise
podo_GLSdata <- snakemake@input$genewise

## --- dGLS

GLSdata  <- read.table(podo_GLSdata, sep="\t", header=FALSE, stringsAsFactors=FALSE)
# Remove redundant columns
GLSdata <- GLSdata[,c(1,2,4,6)]
names(GLSdata) <- c("gene_id", "Tree", "tr1_log.likelihood", "tr2_log.likelihood")
# sapply(GLSdata, class)

# Add difference in likelihood
GLSdata$diff <- (GLSdata$tr1_log.likelihood - GLSdata$tr2_log.likelihood)

# Plot
p_dgls <- ggplot(GLSdata, aes(x=gene_id, y=diff, fill=Tree)) + geom_bar(stat="identity", position="identity") + 
  ylab(expression(Delta~GLS)) + xlab("Marker") +
  theme_bw() +
  theme(legend.position="none") 

## --- SLS per gene

SLSdata_gene  <- read.table(podo_SLSdata_gene, sep="\t", header=FALSE, stringsAsFactors=FALSE)

# Remove redundant columns
SLSdata_gene <- SLSdata_gene[,c(1,3,5)]
names(SLSdata_gene) <- c("gene_id", "AB", "AC") # The clades in the Podosporaceae family T1 = AB, T2 = AC
SLSdata_gene <- gather(SLSdata_gene, "Tree", "Sites", -gene_id)


p_slsg <- ggplot(SLSdata_gene, aes(x=gene_id, y=Sites, fill=Tree)) +
  geom_col(aes(fill = Tree), width = 0.7) + 
  ylab("Sites") + xlab("Marker") +
  theme_bw()


# Put them together
together <- plot_grid(p_dgls, p_slsg, labels = c('A', 'B'), rel_widths = c(0.9, 1.3))
# Save it
ggsave(plot = together, snakemake@output[[1]], width = 6, height = 2.5)

