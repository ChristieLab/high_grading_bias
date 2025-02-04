### Ascertainment Bias in expression data 
### Randomly assign treatments using common garden whelks to test the effects of high gradning bias 
### written by andy lee
### last edited 8/9/2024
library(DESeq2)
library("rtracklayer")
library(ggplot2)

colours <- khroma::color("batlow")
manual_colors <- colours(4, range=c(0.1,0.8))

gene.counts <- read.table("~/KW/common_garden_cross/3_deg_pipeline/rnasapdes_annotated/featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("~/KW/common_garden_cross/4_DESeq2/annotated_rnaspades/kw_sample_info.csv")


rownames(samples) <- samples$sample_ID #CGE
all(rownames(samples) == colnames(gene.counts)) # should be TRUE 


site <- which(samples$site == "MON") # use the largest pop from Lee et al 2024 
site.meta   <- samples[site, ] 
site.genecounts <- gene.counts[, site]

all(rownames(site.meta) %in% colnames(site.genecounts)) # check sample names are in the same order in the gene count sample info table
all(rownames(site.meta) == colnames(site.genecounts)) # must  be TRUE or do not proceed

write.csv(site.meta, "data/cge_mon_sitemeta.csv")
write.table(site.genecounts, "data/cge_mon_genecounts.txt")


read.table("data/cge_mon_genecounts.txt") 
# create a function to run DESeq2 
run_deseq <- function(x){
  set.seed(x)
  site.meta$treatment <- sample(LETTERS[1:2], nrow(site.meta), TRUE) # randomly assign treatments
  dds.site <- DESeqDataSetFromMatrix(countData = site.genecounts,
                                     colData = site.meta,
                                     design = ~treatment) 
  
  dds.site$treatment <- relevel(dds.site$treatment, "A")
  
  # get gene names 
  assembly <- readGFF("~/KW/common_garden_cross/4_DESeq2/annotated_rnaspades/stringtie_all_merged_kw.gtf")             ## read in merged GTF
  gene_idx <- match(dds.site@rowRanges@partitioning@NAMES, assembly$gene_id)
  
  # create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
  contig.name      <- as.character(assembly$seq[gene_idx])
  transcript_start <- assembly$start[gene_idx]
  transcript_end   <- assembly$end[gene_idx]
  transcript.id  <- assembly$transcript_id[gene_idx]
  gene_names     <- cbind(contig.name, transcript_start, transcript_end, transcript.id)
  
  dds.site@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3],gene_names[,4], sep=",") # adds unique gene names to dds
  which(duplicated(dds.site@rowRanges@partitioning@NAMES)) # check to make sure that no gene names in dds are duplicates
  
  ## Run DESeq 
  dds.site <- DESeq(dds.site)
  return(dds.site)
}

# run deseq with 5 different set.seed randomly chosen 

run1 <- run_deseq(0420)

## get "top 1000 loci" 
get_top_loci <- function(dds){
  res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
  top_lfc2 <- dds[order(abs(res$log2FoldChange), decreasing = TRUE), ][1:1000,]
  top_lfc_vst <- varianceStabilizingTransformation(top_lfc2)
  return(top_lfc_vst)
}

### get statistically significant genes 
get_sig_genes <- function(dds){
  res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
  sig_cge_genes <- dds[which(res$padj < 0.05),]
  print(res[which(res$padj < 0.05),])
  cge_sig_vst <- varianceStabilizingTransformation(sig_cge_genes)
}


run1_vst <- varianceStabilizingTransformation(run1)
run1_top1k_vst <- get_top_loci(run1)
run1_sig_vst <- get_sig_genes(run1)


## plot PCAs 
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

p1 <- plot_pcas(run1_vst)
p2 <- plot_pcas(run1_top1k_vst)
p3 <- plot_pcas(run1_sig_vst)

ggsave(
  "testplot.svg", 
  cowplot::plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"), 
    p3, 
    nrow = 1, 
    labels=c("All DEGs","Top 1000 DEGs", "Outlier DEGs"), 
    align="hv"
  )
)
