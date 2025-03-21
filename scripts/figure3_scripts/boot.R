library(ggplot2); library(snpR)
source("functions.R")

args <- commandArgs(TRUE)
d <- as.character(args[1])
outfile <- as.character(args[2])
iter <- as.numeric(args[3])

d <- readRDS(d)

set.seed(64720 + iter) 
cat("iter:", iter, "\t", "n SNPs:", nsnps(d), "\n")

tres <- run_bootstrapping(d, "pop", n = 10, par = 10, store_pca = FALSE, fst_cut = 0.95)

saveRDS(tres, paste0(outfile, "_boot_r_", iter, ".RDS"))
