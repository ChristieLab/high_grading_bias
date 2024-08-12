### Script to test effects of picking top loci using various clustering programs 
### Written by Andy Lee and Will Hemstrom 
### Last edited 8/8/2024
library(snpR)
library(khroma)
library(cowplot)
library(ggplot2)

panmictic <- readRDS("../ascertainment_bias/data/ms_out_1.RDS")  
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
p1 <- plot_clusters(panmictic_sub, "pop", alt.palette = manual_colors, plot_type = "tSNE")
p2 <- plot_clusters(panmictic_sub, "pop", alt.palette = manual_colors, plot_type = "umap")
p3 <- plot_clusters(panmictic_sub, "pop", alt.palette = manual_colors, plot_type = "dapc")
p4 <- plot_structure(panmictic_sub, "pop", alt.palette = manual_colors, method = "snmf", k = 2:4)

### high fst 
p6 <- plot_clusters(panmictic.top, "pop", alt.palette = manual_colors,  plot_type = "tSNE")
p7 <- plot_clusters(panmictic.top, "pop", alt.palette = manual_colors,  plot_type = "umap")
p8 <- plot_clusters(panmictic.top, "pop", alt.palette = manual_colors,  plot_type = "dapc") #300, 4, 300, 3, 
p9 <- plot_structure(panmictic.top, "pop", alt.palette = manual_colors, method = "snmf", k = 2:4)



fig <- plot_grid(p1$plots$tsne + theme(legend.position="none"), 
          p6$plots$tsne + theme(legend.position="none"), 
          p2$plots$umap + theme(legend.position="none"), 
          p7$plots$umap + theme(legend.position="none"), 
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

save_plot("~/downloads/si.pdf", fig,  base_height = 15, base_width = 11)
