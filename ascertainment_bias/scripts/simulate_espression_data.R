## Simulate RNA-Seq data for DESeq2 

## code to simulate RNA Seq data to evaluate highgrading bias 
## written by Andy Lee 8/6/2024, allrights reserved
## contact at andymuanlee@ gmail.com 


# BiocManager::install("PROPER")
library(PROPER)

## set options for simulation 
ngenes        <- 1000
seqDepth      <- 10000
# lBaselineExpr <- 10
lOD           <- "maqc"  
p.DE          <- 0.05 #5% of genes are DE
# lfc           <- fun.lfc 
sim.seed      <- 4

simOptions <- RNAseq.SimOptions.2grp(ngenes, seqDepth, lOD, p.DE, sim.seed)

simRNAseq(simOptions, n1=5, n2=5)
