### Script to test effects of picking top loci using various clustering programs 
### Written by Andy Lee and Will Hemstrom 
### Last edited 8/8/2024
library(snpR)
library(khroma)
library(cowplot)
library(ggplot2)

panmictic <- readRDS("../data/ms_out_1.RDS")  
panmictic <- calc_global_fst(panmictic, "pop")
fst <- get.snpR.stats(panmictic, "pop", "fst")
high_fst <- which(fst$pairwise$fst >= quantile(fst$pairwise$fst, .95, na.rm = TRUE))

set.seed(527615)
samp <- sample(nrow(panmictic), 15000, replace = FALSE)
panmictic_sub <- panmictic[samp, ]

panmictic.top <- panmictic[high_fst,]
saveRDS(panmictic.top <- panmictic[high_fst,], "ms_out_1_highfst.RDS")
# panmictic.top <- calc_global_fst(panmictic.top, "pop")


## Set colors
colours <- khroma::color("batlow")
manual_colors <- colours(4, range=c(0.1,0.8))

### all SNPs 
p1 <- plot_clusters(panmictic_sub, "pop", alt.palette = manual_colors, plot_type = c("tSNE", "umap"))
p3 <- plot_clusters(panmictic_sub, "pop", alt.palette = manual_colors, plot_type = "dapc", dapc_clustering_max_n_clust = NULL,
                    dapc_clustering_npca = 100, dapc_clustering_nclust = 4, dapc_npca = 100, dapc_ndisc = 3) #100, 4, 100, 3
p4 <- plot_structure(panmictic_sub, "pop", alt.palette = manual_colors, method = "snmf", k = 4)

### high fst 
p6 <- plot_clusters(panmictic.top, "pop", alt.palette = manual_colors,  plot_type = c("tSNE", "umap"))
p8 <- plot_clusters(panmictic.top, "pop", alt.palette = manual_colors,  plot_type = "dapc", dapc_clustering_max_n_clust = NULL, 
                    dapc_clustering_npca = 150, dapc_clustering_nclust = 4, dapc_npca = 150, dapc_ndisc = 3) # 50% of var explained
p9 <- plot_structure(panmictic.top, "pop", alt.palette = manual_colors, method = "snmf", k = 4)



fig <- plot_grid(p1$plots$tsne + theme(legend.position="none"), 
          p6$plots$tsne + theme(legend.position="none"), 
          p1$plots$umap + theme(legend.position="none"), 
          p6$plots$umap + theme(legend.position="none"), 
          p3$plots$dapc + theme(legend.position="none"), 
          p8$plots$dapc + theme(legend.position="none"), 
          p4$plot + theme(legend.position="none"),       
          p9$plot + theme(legend.position="none"), 
        
          ncol = 2, 
          align="hv", 
          labels=c("tSNE", "", "umap", "", "DAPC", "", "sNMF", ""),
          label_size = 10, 
          hjust = 0,
          axis="tbl")  

save_plot("../results/Figure_S1.pdf", fig,  base_height = 15, base_width = 11)
