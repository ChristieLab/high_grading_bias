## code associated with Figure 1 of the ascertainment bias manuscript 
## illustrating the problem with picking high FST SNPs using a pamictic monarch dataset and randomly assigned population codes

## Written by William Hemstrom and Andy Lee XX/XX/XXXX

## set wd and load packages 
setwd("~/assignment_tests/ascertainment_bias/data")
library("snpR")
library("ggplot2")
library(tidyverse)


monarch_colors <- c("#d4a85e","#eb5d2e","#829b51","#393e42")

## read in cleaned up RDS objects
monarchs_random_allsnps <- readRDS("monarchs_random.RDS")
monarchs_random_highfst <- readRDS("monarchs_random_highfst.RDS")

## plot PCAs for figure 1 
pca_all  <- plot_clusters(monarchs_random_allsnps, "pop", "pca", alt.palette = monarch_colors)
pca_high <- plot_clusters(monarchs_random_highfst, "pop", "pca", alt.palette = monarch_colors)


## plot STRUCTURE for fig 1 
struc <- plot_structure(monarchs_high, "pop", 
                        method = "structure", 
                        structure_path = "/usr/bin/structure.exe",
                        k = 2:4, iterations = 500, burnin = 100)


## plot RUBIAS results for figure 1 
rubias_all  <- read.table("../results/rubias/monarchs_random_rubias.txt", header = 1)
rubias_high <- read.table("../results/rubias/monarchs_random_highfst_rubias.txt", header = 1)

ggplot(rubias_all) +
  geom_bar(aes(fill=inferred_collection, x=collection), position="fill") + 
  labs(x= "Site of Origin",  y="Proportion of Inferred Collection") +
  scale_fill_manual(values=monarch_colors[3:4]) +
  theme_bw()
  
ggplot(rubias_high) +
  geom_bar(aes(fill=inferred_collection, x=collection), position="fill") + 
  labs(x= "Site of Origin",  y="Proportion of Inferred Collection") +
  scale_fill_manual(values=monarch_colors)+ theme_bw()


## plot other methods for SI 
plot_clusters(monarchs_high, "pop", "tsne")

struc <- plot_structure(monarchs_high, "pop", 
                        method = "snmf", 
                        structure_path = "/usr/bin/structure.exe",
                        k = 2:4, iterations = 500, burnin = 100)

struc <- plot_structure(monarchs_high, "pop", 
                        method = "snapclust", 
                        structure_path = "/usr/bin/structure.exe",
                        k = 2:4, iterations = 500, burnin = 100)
