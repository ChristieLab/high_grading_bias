library(snpR)

monarchs <- readRDS("data/monarch_nomaf.RDS")
monarchs <- monarchs[pop = "NAM"]

meta <- sample.meta(monarchs)
meta$real_pop<- ifelse(meta$samp == "NAM_M9.10" | grepl("_Mexico", meta$samp),
                       "ENA", "WNA")
sample.meta(monarchs)$pop <- meta$real_pop
monarchs <- filter_snps(monarchs, maf = 0.05)
monarchs <- monarchs[sample(nrow(monarchs),1000, FALSE),]

nboots <- 10

generate_summary_stats <- function(tdata, fst_cut = .95){
  tdata <- calc_global_fst(tdata, "pop")
  fst <- get.snpR.stats(tdata, "pop", "fst")
  high_fst <- which(fst$pairwise$fst >= quantile(fst$pairwise$fst, fst_cut, na.rm = TRUE))
  high_fst <- snpR:::.paste.by.facet(fst$pairwise, c("group", "position"), "_")[high_fst]
  high_fst <- match(high_fst,
                    snpR:::.paste.by.facet(snp.meta(tdata), c("group", "position"), "_"))
  
  tdata_high <- tdata[high_fst,]
  tdata_high <- calc_global_fst(tdata_high, "pop")
  snpR:::.make_it_quiet(tpca <- plot_clusters(tdata_high, "pop", simplify_output = TRUE,)$pca$data)
  tfst <- get.snpR.stats(tdata_high, "pop", "fst")$weighted.means$mean_fst
  
  tpca <- as.data.table(tpca)
  
  cols <- c("PC1", "PC2")
  tpca[,c("pop_means_PC1", "pop_means_PC2")  := lapply(.SD, mean), by = pop, .SDcols = cols]
  tpca[,c("global_means_PC1", "global_means_PC2") := lapply(.SD, mean), .SDcols = cols]
  
  tpca[,global_dist := sqrt((PC1 - global_means_PC1)^2  + (PC2 - global_means_PC2)^2)]
  tpca[,pop_dist := sqrt((PC1 - pop_means_PC1)^2  + (PC2 - pop_means_PC2)^2)]
  sum_global <- sum(tpca$global_dist)
  sum_pop <- sum(tpca$pop_dist)
  
  test_statistic <- sum_pop/sum_global
  
  return(list(test_statistic = test_statistic, tfst = tfst))
}

do_boots <- function(x, nboots){
  out <- data.table(boot = 1:nboots)
  out$test_stat <- 0
  out$fst <- 0
  
  for(i in 1:nboots){
    cat(i, "\n")
    tdata <- x
    sample.meta(tdata)$pop <- sample(sample.meta(tdata)$pop, nsamps(tdata), TRUE)
    
    tout <- generate_summary_stats(tdata, .95)
    out[i, test_stat := tout$test_statistic]
    out[i, fst := tout$tfst]
  }
  
  return(out)
}

get_p_values <- function(observed, null, h0 = c("less", "greater")){
  ec_dists <- lapply(null, ecdf)
  
  p <- numeric(length(observed))
  names(p) <- names(observed)
  for(i in 1:length(ec_dists)){
    p[i] <- ec_dists[[i]](observed[[i]])
    if(h0[i] == "greater"){
      p[i] <- 1-p[i]
    }
  }
  return(p)
}


boots_NAM <- do_boots(monarchs, nboots)
real_NAM <- generate_summary_stats(monarchs, .95)
p_NAM <- get_p_values(real_NAM, boots_NAM[,-1])

# try with real pops
gr <- readRDS("data/monarch_nomaf.RDS")
gr <- gr[pop = c("GUA", "ROT")]
gr <- filter_snps(gr, maf = 0.05)
gr <- gr[sample(nrow(gr),1000, FALSE),]

boots <- do_boots(gr, nboots)
real <- generate_summary_stats(gr, .95)
p <- get_p_values(real, boots[,-1])
