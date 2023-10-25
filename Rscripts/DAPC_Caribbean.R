################################################

### Test for genetic structure
### Authors: Angelo Soto-Centeno & Pedro Ivo Monico
# date: 25 of October, 2023

# load packages
library(adegenet)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(igraph)
library(ggplot2)
library(reshape2)
library(tidyr)

################################################ 

# Clear global environment
rm(list = ls())

#Set working directory with the downloaded data 
setwd("/Users/Pedro/Desktop/Eptesicus_fuscus/phylogenetics/DAPC")

#Upload dataset
Lgetu.VCF <- read.vcfR("Lgetu.vcf")
Lgetu.VCF
        # note: vcf files do not contain population information

# Load population data
pop.data <- read.table("/Users/Pedro/Desktop/Eptesicus_fuscus/phylogenetics/DAPC/pop_assign.txt", sep = "\t", header = T)

# check that vcf samples & pop data match
all(colnames(Lgetu.VCF@gt)[-1] == pop.data$AccessID)
# [1] TRUE

# convert vcf data to genlight
gl.Lgetu <- vcfR2genlight(Lgetu.VCF)
  # you may get a Warning of loci with >2 alleles will be ommited. that's OK

# set ploidy number 
ploidy(gl.Lgetu) <- 2

# add population State designations
pop(gl.Lgetu) <- pop.data$locality

# verify the genlight object that contains the filtered VCF
gl.Lgetu

################################################ 

# Distance tree
  # this uses a UPGMA method and 25 bootstrap reps
tree <- aboot(gl.Lgetu, tree = "upgma", distance = bitwise.dist, sample = 25, showtree = F, cutoff = 50, quiet = T)

cols <- brewer.pal(n = nPop(gl.Lgetu), name = "Dark2")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.Lgetu)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8, font = 3, xpd = TRUE)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")


################################################ 

# Minimum Spanning Network #
  # create clusters of multilocus genotypes (MLG) by genetic distance

Lgetu.dist <- bitwise.dist(gl.Lgetu)
Lgetu.msn <- poppr.msn(gl.Lgetu, Lgetu.dist, showplot = FALSE, include.ties = T)

node.size <- rep(2, times = nInd(gl.Lgetu))
names(node.size) <- indNames(gl.Lgetu)
vertex.attributes(Lgetu.msn$graph)$size <- node.size

set.seed(5)

plot_poppr_msn(gl.Lgetu, Lgetu.msn , palette = brewer.pal(n = nPop(gl.Lgetu), name = "Dark2"), gadj = 70)

################################################ 

## Principal Component Analysis ##
  # convert SNP data into linearly uncorrelated Principal Components to summarize variation between samples

Lgetu.pca <- glPca(gl.Lgetu, nf = 2)
# screeplot of eigenvalues
barplot(100*Lgetu.pca$eig/sum(Lgetu.pca$eig), col = heat.colors(50), main = "PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

# scatterplot of PCs
# convert scores to a dataframe that works with ggplot
Lgetu.pca.scores <- as.data.frame(Lgetu.pca$scores)
Lgetu.pca.scores$pop <- pop(gl.Lgetu)

set.seed(9)

# step by step plot
p <- ggplot(Lgetu.pca.scores, aes(x = PC1, y = PC2, colour = pop))
p <- p + geom_point(size = 2)
p <- p + stat_ellipse(level = 0.95, size = 1) 
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
# print it
p

################################################ 

## DAPC ##
  # maximize the variance among populations in the sample 

pnw.dapc <- dapc(gl.Lgetu, n.pca = 2, n.da = 2)

# scatterplot of DAPC
scatter(pnw.dapc, cex = 2, legend = TRUE, clabel = F, posi.leg = "topright", 
        scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)

## probability of population membership ##
# compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(pnw.dapc, posi = 'top')

# convert the DAPC data into a data.frame
dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.Lgetu)
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

# check data.frame
head(dapc.results, n = 6)

# rename columns
colnames(dapc.results) <- c("Original_Pop", "Sample", "Assigned_Pop", "Posterior_membership_probability")

# plot dpac.results of the reorganized pivot_longer data.frame into population membership
p <- ggplot(dapc.results, aes(x = Sample, y = Posterior_membership_probability, fill = Assigned_Pop))
p <- p + geom_bar(stat = 'identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

################################################ 
#END 