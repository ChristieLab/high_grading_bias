### Script to test effects of picking top loci using various clustering programs 
### Written by Andy Lee and Will Hemstrom 
### Last edited 8/8/2024
library(snpR)
library(khroma)
library(cowplot)

panmictic <- readRDS("../ascertainment_bias/data/ms_out_1.RDS")  
panmictic <- calc_global_fst(panmictic, "pop")
fst <- get.snpR.stats(panmictic, "pop", "fst")
high_fst <- which(fst$pairwise$fst >= quantile(fst$pairwise$fst, .95, na.rm = TRUE))

panmictic.top <- panmictic[high_fst,]
saveRDS(panmictic.top <- panmictic[high_fst,], "ms_out_1_highfst.RDS")
# panmictic.top <- calc_global_fst(panmictic.top, "pop")


## Set colors
colours <- color("batlow")
manual_colors <- colours(4, range=c(0.1,0.8))

### all SNPs 
p1 <- plot_clusters(panmictic, "pop", alt.palette = manual_colors, method = "tSNE")
p2 <- plot_clusters(panmictic, "pop", alt.palette = manual_colors, method = "umap")
p3 <- plot_clusters(panmictic, "pop", alt.palette = manual_colors, method = "dapc")
p4 <- plot_structure(panmictic, "pop", alt.palette = manual_colors, method = "snmf", k = 2:4, iterations =100000)

### high fst 
p6 <- plot_clusters(panmictic.top, "pop", alt.palette = manual_colors,  plot_type = "tSNE")
p7 <- plot_clusters(panmictic.top, "pop", alt.palette = manual_colors,  plot_type = "umap")
p8 <- plot_clusters(panmictic.top, "pop", alt.palette = manual_colors,  plot_type = "dapc")
p9 <- plot_structure(panmictic.top, "pop", alt.palette = manual_colors, method = "snmf", k = 2:4, iterations = 100000)


plot_grid(#p1$plots$tsne, 
          p6$plots$tsne, 
          #p2$plots$umap, 
          p7$plots$umap, 
          #p3$plots$dapc, 
          p8$plots$dapc, 
          #p4$plot,       
          p9$plot, 
        
          ncol = 1, 
          align="hv", 
          #labels=c("tSNE", "umap", "DAPC", "sNMF", "", "", "", ""),
          # vjust = -0.5, 
          axis="tbl")  
