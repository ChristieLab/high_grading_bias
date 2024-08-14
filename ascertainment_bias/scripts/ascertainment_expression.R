### Ascertainment Bias in expression data 
### Randomly assign treatments using common garden whelks to test the effects of high gradning bias 
### written by andy lee
### last edited 8/9/2024
library(DESeq2)
library("rtracklayer")
library(ggplot2)


#gene.counts <- read.table("~/KW/common_garden_cross/3_deg_pipeline/rnasapdes_annotated/featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
gene.counts <- read.table("~/tae/data/deseq2/tae_genome_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)

drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

#samples <- read.csv("~/KW/common_garden_cross/4_DESeq2/annotated_rnaspades/kw_sample_info.csv")
samples <- read.csv("~/tae/data/tae_sample_info.csv")

#rownames(samples) <- samples$sample_ID
rownames(samples) <- samples$sample_bam
all(rownames(samples) == colnames(gene.counts)) # should be TRUE 

set.seed(981934)
samples$treatment <- sample(LETTERS[1:2], nrow(samples), TRUE) # randomly assign treatments


### examine response to temperature in each site -----------------------------##
#----Sample processing ----#

site <- which(samples$site == "MON") # use the largest pop from Lee et al 2024 

# site <- sample(site, 6)
site.meta   <- samples[site, ] 
site.genecounts <- gene.counts[, site]

all(rownames(site.meta) %in% colnames(site.genecounts)) # check sample names are in the same order in the gene count sample info table
all(rownames(site.meta) == colnames(site.genecounts)) # must  be TRUE or do not proceed

dds.site <- DESeqDataSetFromMatrix(countData = site.genecounts,
                                   colData = site.meta,
                                   design = ~treatment) 

dds.site$treatment <- relevel(dds.site$treatment, "A")
dds.site

## gene names 
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)

#assembly <- readGFF("~/KW/common_garden_cross/4_DESeq2/annotated_rnaspades/stringtie_all_merged_kw.gtf")             ## read in merged GTF
assembly <- readGFF("~/tae/data/tae_all_merged_genome.gtf")             ## read in merged GTF

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
res.site <- results(dds.site, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res.site[, 6] < 0.05))
summary(res.site)

## TEST CODE
# top_lfc  <- which(abs(res.site$log2FoldChange) > quantile(abs(res.site$log2FoldChange), .999, na.rm = T)) 

top_lfc2 <- dds.site[order(abs(res.site$log2FoldChange), decreasing = TRUE), ][1:1000,]
top_lfc_vst <- varianceStabilizingTransformation(top_lfc2)

plot.all.data <- plotPCA(top_lfc_vst, intgroup = "treatment", returnData=TRUE)
percentVar <- round(100*attr(plot.all.data, "percentVar"))

## color PCA by treatment
p_cge <- ggplot(plot.all.data, aes(PC1, PC2, color=treatment, shape=treatment)) +
  # ggtitle(paste("Top",num.genes,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +
  scale_color_manual(values=c("royalblue3","tomato1")) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

cge_sig_vst <- varianceStabilizingTransformation(dds.site)

### plot with only significant DEGs 
sig_cge_genes <- dds.site[which(res.site$padj < 0.05),]
cge_sig_vst <- varianceStabilizingTransformation(sig_cge_genes)
plot.cge_sig <- plotPCA(cge_sig_vst, intgroup = "treatment", returnData=TRUE)
percentVar <- round(100*attr(plot.cge_sig, "percentVar"))

p_sig_cge <- ggplot(plot.cge_sig, aes(PC1, PC2, color=treatment, shape=treatment)) +
  # ggtitle(paste("Top",num.genes,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +
  scale_color_manual(values=c("royalblue3","tomato1")) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()


## tae/mon sig genes 
sig_cge_genes <- dds.site[which(res.site$padj < 0.05),]
cge_sig_vst <- varianceStabilizingTransformation(sig_cge_genes)
plot.cge_sig <- plotPCA(cge_sig_vst, intgroup = "treatment", returnData=TRUE)
percentVar <- round(100*attr(plot.cge_sig, "percentVar"))

p_sig_cge <- ggplot(plot.cge_sig, aes(PC1, PC2, color=treatment, shape=treatment)) +
  # ggtitle(paste("Top",num.genes,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +
  scale_color_manual(values=c("royalblue3","tomato1")) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()


par(mar = c(0,0,0,0))


cowplot::plot_grid(
  p_cge + theme(legend.position="none"),
  p_tae + 
    labs(color="Randomly assigned treatment", shape="Randomly assigned treatment")
    + theme(legend.position="none"),  
  p_sig_cge + theme(legend.position="none"),
  ncol = 2, 
  # labels=c("No existing treatments", "Randomly Shuffled treatments", "Significant DEGs"), 
  align="h"
  )


# plot all 4 tae sites
# cowplot::plot_grid(
#   p_mon + theme(legend.position="none"),
#   p_dic + theme(legend.position="none"),
#   p_nap + theme(legend.position="none"),
#   p_pol + theme(legend.position="none"), 
#   ncol = 2, 
#   align="hv", 
#   labels=c("mon", "dic", "nap", "pol"),
#   # label_size = 10, 
#   # hjust = 0,
#   axis="trbl") 

