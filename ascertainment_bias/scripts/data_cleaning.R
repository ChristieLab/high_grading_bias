### Script written by Andy Lee 
### Last Modified: April 2024
### Standardizing datasets for MER population assignment power analysis proposal 
### using snpR package 

### Load packages and setwd 
library("snpR")
setwd("~/assignment_tests/ascertainment_bias/data")

## Standardize each data set 



###############################
#### Testing high FST SNPs ####
###############################


#### Monarch data for testing ascertainment bias using highly and non differentiated populations
get_high_fst <- function(rds, percent_highest){
  rds <- calc_global_fst(rds, "pop")
  fst <- get.snpR.stats(rds, "pop", "fst")
  high_fst <- which(fst$pairwise$fst >= quantile(fst$pairwise$fst, percent_highest, na.rm = TRUE))
  high_fst <- snpR:::.paste.by.facet(fst$pairwise, c("group", "position"), "_")[high_fst]
  high_fst <- match(high_fst,
                    snpR:::.paste.by.facet(snp.meta(rds), c("group", "position"), "_"))
  high_fst_rds <- rds[high_fst,]
  return(high_fst_rds)
}

monarchs <- readRDS("../../assignment_tests/assignment_tests/data/monarch_nomaf.RDS")
colnames(sample.meta(monarchs)) <- c("sampID", "pop", ".sample.id")
sample.meta(monarchs)$sampID <- gsub("\\.", "_", sample.meta(monarchs)$sampID)
sample.meta(monarchs) <- sample.meta(monarchs)[,1:2]

## panmictic north america samples 
monarch_nam <- filter_snps(monarchs[pop = "NAM"], maf = 0.05, maf_facets= "pop")
meta <- sample.meta(monarch_nam)
sample.meta(monarch_nam)$pop <- ifelse(meta$samp == "NAM_M9.10" | grepl("_Mexico", meta$samp), "ENA", "WNA")
monarch_nam_highfst <- get_high_fst(monarch_nam, .95)
saveRDS(monarch_nam_highfst, "monarch_nam_highfst.RDS")

## highly differentiated GUA and ROT samples 
monarch_gr <- filter_snps(monarchs[pop = c("GUA", "ROT")], maf = 0.05, maf_facets= "pop")
monarch_gr_highfst <- get_high_fst(monarch_gr, .95)
saveRDS(monarch_gr, "monarch_gr_highfst.RDS")
