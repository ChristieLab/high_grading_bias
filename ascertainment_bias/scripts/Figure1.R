## code associated with Figure 1 of the high grading bias manuscript 
## illustrating the problem with picking high FST SNPs using a panmictic monarch dataset and randomly assigned population codes

## Written by William Hemstrom and Andy Lee XX/XX/XXXX

## set wd and load packages 
setwd("~/assignment_tests/ascertainment_bias/data")
library("snpR")
library("ggplot2")
library(tidyverse)

## read in cleaned up RDS objects
monarchs_random_allsnps <- readRDS("monarchs_random.RDS")
monarchs_random_highfst <- readRDS("monarchs_random_highfst.RDS")

## Set colors
monarch_colors <- c("#d4a85e","#eb5d2e","#829b51","#393e42")

## plot PCAs for figure 1 
pca_all  <- plot_clusters(monarchs_random_allsnps, "pop", "pca", alt.palette = monarch_colors)

pca_all$plots$pca <- pca_all$plots$pca + 
  guides(color=guide_legend(title = stringr::str_wrap("Randomly Assigned Population", width =20)))

pca_high <- plot_clusters(monarchs_random_highfst, "pop", "pca", alt.palette = monarch_colors)


## plot STRUCTURE for fig 1 
setwd("~/assignment_tests/ascertainment_bias/results/structure/all_sites_structure_res/")
structure_allsnps <- plot_structure("merged", k = 4, facet = sample.meta(monarchs_random_allsnps)$pop, alt.palette = monarch_colors, facet.order = c("A", "B", "C", "D"))

setwd("~/assignment_tests/ascertainment_bias/results/structure/high_fst_structure_res/")
structure_highfst <- plot_structure("merged", k = 4, facet = sample.meta(monarchs_random_highfst)$pop, alt.palette = monarch_colors, facet.order = c("A", "B", "C", "D"))

## plot RUBIAS results for figure 1 
rubias_all  <- read.table("~/assignment_tests/ascertainment_bias/results/rubias/monarchs_random_rubias.txt", header = 1)
rubias_high <- read.table("~/assignment_tests/ascertainment_bias/results/rubias/monarchs_random_highfst_rubias.txt", header = 1)

p_rubias_all <- ggplot(rubias_all) +
  geom_bar(aes(fill=inferred_collection, x=collection), position="fill", width =0.975) + 
  labs(x= "Randomly Assigned Population",  y="Proportion of Individuals Assigned") +
  scale_fill_manual(values=monarch_colors[3:4]) +
  theme_bw() + 
  guides(fill=guide_legend(title = stringr::str_wrap("Inferred Population", width =20)))
  
p_rubias_high <- ggplot(rubias_high) +
  geom_bar(aes(fill=inferred_collection, x=collection), position="fill", width =0.975) + 
  labs(x= "Randomly Assigned Population",  y="Proportion of Individuals Assigned") +
  scale_fill_manual(values=monarch_colors)+ theme_bw() + 
  guides(fill=guide_legend(title = stringr::str_wrap("Inferred Population", width =20)))

library(cowplot)
library(ggpubr)
top_lg <- get_legend(pca_all$plots$pca + 
                       theme(legend.title = element_text(size=14),
                             legend.text = element_text(size=12)))
mid_lg <- get_legend(structure_allsnps$plot + 
                       theme(legend.title = element_text(size=14),
                             legend.text = element_text(size=12)))
btm_lg <- get_legend(p_rubias_high + 
                       theme(legend.title = element_text(size=14),
                             legend.text = element_text(size=12)))

figure_1 <- plot_grid(
  pca_all$plots$pca + 
    guides(color="none") + 
    ggtitle("All SNPs") + 
    theme(plot.title = element_text(size=20), 
          axis.text = element_text(angle=0, size=12), 
          axis.title = element_text(size=14)), 
  
  pca_high$plots$pca + guides(color="none") + 
    ggtitle(expression('Top' * ' ' * '5%' * ' ' *italic("F")[ST])) + 
    theme(plot.title = element_text(size=20), 
          axis.text = element_text(angle=0, size=12), 
          axis.title = element_text(size=14)), 
  top_lg,
  
  structure_allsnps$plot + 
    guides(fill="none", color="none") + 
    xlab("Randomly Assigned Population") + 
    theme(axis.text = element_text(angle=0, size=12), 
          axis.title = element_text(size=14)),
  
  structure_highfst$plot + 
    guides(fill="none", color="none") + 
    xlab("Randomly Assigned Population") + 
    theme(axis.text = element_text(angle=0, size=12), 
          axis.title = element_text(size=14)),
  mid_lg,
  
  p_rubias_all + guides(fill="none") + 
    scale_x_discrete(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    theme(axis.text = element_text(angle=0, size=12), 
          axis.title = element_text(size=14)),
  
  p_rubias_high + guides(fill="none") + 
    scale_x_discrete(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    theme(axis.text = element_text(angle=0, size=12),
          axis.title = element_text(size=14)), 
  
  btm_lg,
  ncol = 3, axis = "tr", align="hv", rel_widths =c(1,1,.3)
)

save_plot("~/assignment_tests/ascertainment_bias/results/Figure1.PDF", figure_1,  base_height = 11, base_width = 15)

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
