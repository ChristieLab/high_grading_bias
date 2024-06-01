setwd("~/assignment_tests/ascertainment_bias/data")
library("snpR")
library("ggplot2")

## read in cleaned up RDS objects
monarchs_random_allsnps <- readRDS("monarchs_random.RDS")
monarchs_random_highfst <- readRDS("monarchs_random_highfst.RDS")

## plot PCAs for figure 1 
pca_all  <- plot_clusters(monarchs_random_allsnps, "pop", "pca")
pca_high <- plot_clusters(monarchs_random_highfst, "pop", "pca")


## plot STRUCTURE for fig 1 
struc <- plot_structure(monarchs_high, "pop", 
                        method = "structure", 
                        structure_path = "/usr/bin/structure.exe",
                        k = 2:4, iterations = 500, burnin = 100)

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
