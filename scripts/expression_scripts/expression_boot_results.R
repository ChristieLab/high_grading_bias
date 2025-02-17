library(ggplot2)
library(ggridges)
library(DESeq2)

setwd("~/high_grading_bias/results/expression_hgb/")
temp  <- list.files(pattern = "\\.RDS")
res <-  lapply(temp[1:100], readRDS)

## get f stat f
fstat_top <- res |> purrr::map(1) 
fstat_top <- unlist(fstat_top)

fstat_sig <- res |> purrr::map(2) 
fstat_sig[sapply(fstat_sig, is.null)] <- NA # keeps all the NULL results
fstat_sig <- unlist(fstat_sig)

sum(is.na(fstat_sig) == TRUE) ## How many runs could we not calculate F-stat 
fstat_sig[sapply(fstat_sig, is.na)] <- 0

d1 <- data.frame(fstat = fstat_top, loci = "Top 1000 DEGs")
d2 <- data.frame(fstat = fstat_sig, loci = "Significant DEGs")

df <- rbind(d1,d2)

p1 <- ggplot() + 
  geom_density(data=d1, aes(x=fstat)) + 
  geom_segment(aes(x=65.01611, xend=65.01611, y=0, yend=0.012), color="red") +
  theme_bw() + 
  scale_x_continuous(limits = c(-10,350)) +
  theme(axis.title.y=element_blank())


  
p2 <- ggplot() +
  geom_density(data=d2, aes(x=fstat)) + 
  geom_segment(aes(x=1.067702, xend=1.067702, y=0, yend=0.28), color="red") +
  scale_x_continuous(limits = c(-10,350))  +
  theme_bw() +
  theme(axis.title.y=element_blank())

p1
p2
# ggplot(data=df, aes(x=fstat, y=loci)) +
#   geom_density_ridges(scale = .9, alpha = 0.8) +
#   labs(x = expression(~italic(F)~-statistic), y ="Selected Loci") +
#   scale_y_discrete(expand = c(0, 0)) +
#   coord_cartesian(clip="off") +
#   theme_ridges(center_axis_labels = T) +
#   theme(legend.position = "none", text=element_text(size=20))+
#   theme(axis.title.y=element_blank())
  
## how many runs had significant DEGs? 
setwd("/Users/andy/Library/Application Support/Mountain Duck/Volumes.noindex/Negishi Scratch.localized/ascertainment_bias/results/")
temp  <- list.files(pattern = "\\.txt")
sig_genes <-  lapply(temp[1:100], read.table)

readRDS("deseq_permute_sig_genes_.RDS")

num_outlier <-read.table("/Users/andy/Library/Application Support/Mountain Duck/Volumes.noindex/Negishi Scratch.localized/ascertainment_bias/results/deseq_permute_sig_gene_count.txt")

mean(num_outlier[,2])
median(num_outlier[,2])
max(num_outlier[,2])
plot(density(num_outlier[,2]))
sd(num_outlier[,2])


dplyr::arrange(as.data.frame(num_outlier), V2) ## 5 HAD more than 1000 outliers??
deseq_outliers <- deseq_dims[deseq_dims[,1]<=1000, ] ## remove top outliers 

sum(num_outlier[,2] %in% c(0,1,2)) #how many only ID'd 1-2 sig DEGs
sum(num_outlier[,2] %in% c(0))
sum(num_outlier[,2] %in% c(1))
sum(num_outlier[,2] %in% c(2))
table(num_outlier)

##### plot pcas of the "median" case 
fstat_sig[sapply(fstat_sig, is.na)] <- 0

View(fstat_sig)
mean(fstat_top)
mean(fstat_sig)

d1 <- data.frame(fstat = fstat_top, loci = "Top 1000 DEGs")
d2 <- data.frame(fstat = fstat_sig, loci = "Significant DEGs")

## run 69 had the most representative f-stats
d1[69, ]
d2[69, ]
df <- cbind(d1,d2)
View(df)
## plot example PCA
treatment <- c(rep('A', nrow(site.meta)/2), rep('B', nrow(site.meta)/2))


iter <- 70 ## the 70th run is 69 in the above df
set.seed(930 + iter)

dds <- run_deseq(treatment)

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
dds[which(abs(res$log2FoldChange) > 2),]

sig_cge_genes <- dds[which(res$padj < 0.05),]
print(res[which(res$padj < 0.05),])





top_genes <- get_top_loci(dds)
sig_genes <- get_sig_genes(dds)



fstat_pcas(top_genes)
fstat_pcas(sig_genes)

p3 <- plot_pcas(top_genes[[2]])
p4 <- plot_pcas(sig_genes[[2]])

plot <- cowplot::plot_grid(p1, p3, p2, p4, labels="AUTO", nrow = 2, align = "hv")
ggsave(plot, filename = "../../figures/SI_express_boot.svg", dpi = 400, width = 8, height = 6, units = "in")


## plot them all together 

### functions ====
run_deseq <- function(treatment){
  site.meta$treatment <- sample(treatment, length(treatment), FALSE) # randomly assign treatments
  dds.site <- DESeqDataSetFromMatrix(countData = site.genecounts,
                                     colData = site.meta,
                                     design = ~treatment) 
  
  dds.site$treatment <- relevel(dds.site$treatment, "A")
  
  # # get gene names 
  # assembly <- readGFF("~/KW/common_garden_cross/4_DESeq2/annotated_rnaspades/stringtie_all_merged_kw.gtf")             ## read in merged GTF
  # gene_idx <- match(dds.site@rowRanges@partitioning@NAMES, assembly$gene_id)
  # 
  # # create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
  # contig.name      <- as.character(assembly$seq[gene_idx])
  # transcript_start <- assembly$start[gene_idx]
  # transcript_end   <- assembly$end[gene_idx]
  # transcript.id  <- assembly$transcript_id[gene_idx]
  # gene_names     <- cbind(contig.name, transcript_start, transcript_end, transcript.id)
  # 
  # dds.site@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3],gene_names[,4], sep=",") # adds unique gene names to dds
  # which(duplicated(dds.site@rowRanges@partitioning@NAMES)) # check to make sure that no gene names in dds are duplicates  # which(duplicated(dds.site@rowRanges@partitioning@NAMES)) # check to make sure that no gene names in dds are duplicates
  
  ## Run DESeq 
  dds.site <- DESeq(dds.site)
  return(dds.site)
}
get_top_loci <- function(dds){
  res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
  top_lfc2 <- dds[order(abs(res$log2FoldChange), decreasing = TRUE), ][1:1000,]
  top_lfc_vst <- varianceStabilizingTransformation(top_lfc2)
  return(list(top_lfc2, top_lfc_vst))
} ## get "top 1000 loci" 
get_sig_genes <- function(dds){
  res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
  sig_cge_genes <- dds[which(res$padj < 0.05),]
  
  if(nrow(sig_cge_genes) > 1){
    tryCatch(
      {
        return(list(sig_cge_genes, varianceStabilizingTransformation(sig_cge_genes)))
      },
      error = function(cond) {
        message("Here's the original error message:")
        message(conditionMessage(cond))
        # Choose a return value in case of error
        list(sig_cge_genes, NULL)
      }
    )
  }
  else{
    return(list(sig_cge_genes, NULL))
  }
}
plot_pcas <- function(top_lfc_vst){
  plot.all.data <- plotPCA(top_lfc_vst, intgroup = "treatment", returnData=TRUE)
  percentVar <- round(100*attr(plot.all.data, "percentVar"))
  
  ## plot PCA by treatment using top 1000 DEGs 
  ggplot(plot.all.data, aes(PC1, PC2, color=treatment, shape=treatment)) +
    # ggtitle(paste("Top",num.genes,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +
    scale_color_manual(values=c("#0C335E","#C19039")) + 
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()
}
