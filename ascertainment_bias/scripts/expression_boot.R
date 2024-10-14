### Ascertainment Bias in expression data 
### bootstrap 100x 

### Randomly assign treatments using common garden whelks to test the effects of high gradning bias 
### written by andy lee AND william hemstrom 
### last edited 10/02/2024


## Load packages 

library(DESeq2)
library("rtracklayer")
library(ggplot2)

args <- commandArgs(TRUE)
iter <- as.numeric(args[1])

site.meta <- read.csv("/scratch/negishi/lee3617/ascertainment_bias/data/cge_mon_sitemeta.csv", row.names = 1)
site.genecounts <- read.table("/scratch/negishi/lee3617/ascertainment_bias/data/cge_mon_genecounts.txt") 

# site.meta <- read.csv("../data/cge_mon_sitemeta.csv", row.names = 1)
# site.genecounts <- read.table("../data/cge_mon_genecounts.txt") 

set.seed(930 + iter)
paste0("seed ", 930+iter)

treatment <- sample(LETTERS[1:2], 30, TRUE)
  
# create a function to run DESeq2 

run_deseq <- function(treatment){
  site.meta$treatment <- sample(treatment, length(treatment), TRUE) # randomly assign treatments
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

# get_sig_genes <- function(dds){
#   res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
#   sig_cge_genes <- dds[which(res$padj < 0.05),]
#   
#   
#   if(nrow(sig_cge_genes) > 1){
#     return(list(sig_cge_genes, varianceStabilizingTransformation(sig_cge_genes)))
#   } 
#   else{
#     return(list(sig_cge_genes, NULL))
#   }
# } ### get statistically significant genes 

# 
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

## plot PCAs 
fstat_pcas <- function(vst){
  if(!is.null(vst[[2]])){
    pca.dat <- plotPCA(vst[[2]], intgroup = "treatment", returnData=TRUE)
    Fstat <- summary(manova(as.matrix(pca.dat[,c(1:2)]) ~pca.dat$treatment))$stats[1,"approx F"]
    return(Fstat)
  }
  else{
    return(NULL)
  }
  
  ## plot PCA by treatment using top 1000 DEGs 
  # percentVar <- round(100*attr(plot.all.data, "percentVar"))
  # ggplot(plot.all.data, aes(PC1, PC2, color=treatment, shape=treatment)) +
  #   # ggtitle(paste("Top",num.genes,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +
  #   scale_color_manual(values=c("#0C335E","#C19039")) +
  #   geom_point(size=3) +
  #   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #   coord_fixed()
}

dds <- run_deseq(treatment)
top_genes <- get_top_loci(dds)
sig_genes <- get_sig_genes(dds)

f_top <- fstat_pcas(top_genes)
f_sig <- fstat_pcas(sig_genes)

saveRDS(list(f_top, f_sig, sig_genes[[1]]), paste0("/scratch/negishi/lee3617/ascertainment_bias/results/deseq_boot_res_", iter, ".RDS"))  