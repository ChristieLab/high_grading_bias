### Script written by Andy Lee 
### Last Modified: April 2024
### Standardizing datasets for MER population assignment power analysis proposal 
### using snpR package 

### Load packages and setwd 
# remotes::install_github("hemstrow/snpR", ref = "dev")
library("snpR")
setwd("~/assignment_tests/ascertainment_bias/data")

## Standardize each data set 


#### function to get highest fst snps 
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

#### Monarch data for testing ascertainment bias using highly and non differentiated populations
monarchs_random <- readRDS("monarchs_random.RDS")
monarchs_gr     <- readRDS("monarchs_gr.RDS")

sim_cryptic     <- readRDS("cryptic_genotypes_raw.RDS")
  sample.meta(sim_cryptic)$sampID <- rownames(sample.meta(sim_cryptic)) # standardize sample meta names 
  sample.meta(sim_cryptic) <- sample.meta(sim_cryptic)[,c(3,1)]
  sample.meta(sim_cryptic)$pop <- as.character(sample.meta(sim_cryptic)$pop)
  saveRDS(sim_cryptic,    "sim_cryptic.RDS")
  
monarchs_random_highfst <- get_high_fst(monarchs_random, .95)
monarchs_gr_highfst     <- get_high_fst(monarchs_gr, .95)
sim_cryptic_highfst     <- get_high_fst(sim_cryptic, .95)

saveRDS(monarchs_random_highfst, "monarchs_random_highfst.RDS")
saveRDS(monarchs_gr_highfst,     "monarchs_gr_highfst.RDS")
saveRDS(sim_cryptic_highfst ,    "sim_cryptic_highfst.RDS")

# ## monarchs clean up
# monarchs <- readRDS("../../assignment_tests/data/monarch_nomaf.RDS")
# colnames(sample.meta(monarchs)) <- c("sampID", "pop", ".sample.id")
# sample.meta(monarchs)$sampID <- gsub("\\.", "_", sample.meta(monarchs)$sampID)
# sample.meta(monarchs) <- sample.meta(monarchs)[,1:2]
# 
# ## panmictic north america samples 
# monarch_nam <- filter_snps(monarchs[pop = "NAM"], maf = 0.05, maf_facets= "pop")
# meta <- sample.meta(monarch_nam)
# sample.meta(monarch_nam)$pop <- ifelse(meta$samp == "NAM_M9.10" | grepl("_Mexico", meta$samp), "ENA", # "WNA")
# monarch_nam_highfst <- get_high_fst(monarch_nam, .95)
# saveRDS(monarch_nam_highfst, "monarch_nam_highfst.RDS")
# 
# ## highly differentiated GUA and ROT samples 
# monarch_gr <- filter_snps(monarchs[pop = c("GUA", "ROT")], maf = 0.05, maf_facets= "pop")
# monarch_gr_highfst <- get_high_fst(monarch_gr, .95)
# saveRDS(monarch_gr, "monarch_gr_highfst.RDS")
# 
# 
# ## randomly assign A-D to panmictic north america samples 
# monarchs_random <- filter_snps(monarchs[pop = "NAM"], maf = 0.05)
# sample.meta(monarchs_random)$pop <- sample(c("A", "B", "C", "D"), nsamps(monarchs_random), TRUE)
# saveRDS(monarchs_random, "monarch_nam_random.RDS")