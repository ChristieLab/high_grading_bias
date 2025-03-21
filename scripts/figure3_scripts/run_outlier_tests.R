library(snpR); library(OutFLANK); library(pcadapt)

args <- commandArgs(TRUE)
dat <- readRDS(as.character(args[1]))
outfile <- as.character(args[2])

run_outflank <- function(dat){
  sn <- format_snps(dat, "sn", interpolate = FALSE)
  header_cols <- ncol(snp.meta(dat)) - 1
  sn <- t(sn[,-c(1:header_cols)])
  sn[is.na(sn)] <- 9
  
  pops <- sample.meta(dat)$pop
  
  cat("Preparing data...\n")
  FstDataFrame <- OutFLANK::MakeDiploidFSTMat(sn,
                                              locusNames = snp.meta(dat)$.snp.id,
                                              popNames = pops)
  
  
  cat("Running OutFLANK.\n")
  out1 <- OutFLANK(FstDataFrame, NumberOfSamples = length(unique(pops)))
  out1 <- cbind(out1, snp.meta(dat))
  
  return(out1)
}

run_PCAdapt <- function(dat, K = 2, min.maf = 0.001, LDsize = 200, LDthresh = 0.1){
  tmp <- tempfile(tmpdir = ".")

  if("chr" %in% colnames(snp.meta(dat))){
    format_snps(dat, "plink", "pop", chr = "chr", outfile = tmp)
  }
  else{
    format_snps(dat, "plink", "pop", chr = "group", outfile = tmp)
  }

  pcadpt_input <- read.pcadapt(paste0(tmp, ".bed"), "bed")
  
  pcapt_res <- pcadapt::pcadapt(pcadpt_input, K = K, min.maf = min.maf, LD.clumping = list(size = LDsize, thr = LDthresh))
  
  rmf <- list.files(".", basename(tmp))
  file.remove(rmf)  

  return(cbind(snp.meta(dat), p = pcapt_res$pvalues))
}


dat <- filter_snps(dat, maf = 0.05, maf_facets = "pop")

outflank_res <- run_outflank(dat)
cat("Finished with outFLANK, saving results to:", paste0(outfile, "_outflank_res.RDS"), "\n")
saveRDS(outflank_res, paste0(outfile, "_outflank_res.RDS"))

PCAdapt_res <- run_PCAdapt(dat)
cat("Finished with pcadapt, saving results to:", paste0(outfile, "_PCAdapt_res.RDS"), "\n")
saveRDS(PCAdapt_res, paste0(outfile, "_PCAdapt_res.RDS"))
