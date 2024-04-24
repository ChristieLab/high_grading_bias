### Ascertainment Bias in expression data 
### Randomly shuffle treatments in TAE to test the effects

gene.counts <- read.table("tae_genome_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("../tae_sample_info.csv")
rownames(samples) <- samples$sample_bam
all(rownames(samples) == colnames(gene.counts)) # should be TRUE 

samples$true_treat <- samples$treatment
samples$treatment <- sample(LETTERS[1:2], nrow(samples), TRUE)



### examine response to temperature in each site -----------------------------##
#----Sample processing ----#

site <- which(samples$site == "NAP") 
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
assembly <- readGFF("tae_all_merged_genome.gtf")             ## read in merged GTF
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
top_lfc <- which(abs(res.site$log2FoldChange) > quantile(abs(res.site$log2FoldChange), .95, na.rm = T)) 

top_lfc <- which(abs(res.site$log2FoldChange) > quantile(abs(res.site$log2FoldChange), .99999, na.rm = T)) 

top_lfc_vst <- varianceStabilizingTransformation(dds.site[top_lfc, ])


plot.all.data <- plotPCA(top_lfc_vst, intgroup = "treatment", returnData=TRUE)
percentVar <- round(100*attr(plot.all.data, "percentVar"))

## color PCA by treatment
ggplot(plot.all.data, aes(PC1, PC2, color=treatment, shape=treatment)) +
  # ggtitle(paste("Top",num.genes,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +
  scale_color_manual(values=c("darkblue","firebrick")) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

