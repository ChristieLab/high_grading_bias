### Script written by Andy Lee 
### Last Modified: April 2024
### Standardizing datasets for MER population assignment power analysis proposal 
### using snpR package 

### Load packages and setwd 
library("snpR")
setwd("~/assignment_power_analyses/data/")

## Standardize each data set 

## Kellet's whelks RNA-seq from Lee et al 2024
kw.dat <- readRDS("original_data/kw_rnaseq_allsnps.RDS")
summarize_facets(kw.dat, "pop") # check to see if each pop has more than 5 indivuals 

sample.meta(kw.dat)$sample_ID <- gsub(".sorted.bam", "", sample.meta(kw.dat)$sample_ID) # remove periods 
sample.meta(kw.dat) <-  sample.meta(kw.dat)[,c(1,3)] # retain only necessary meta data 
colnames(sample.meta(kw.dat)) <- c("sampID", "pop", ".sample.id") # rename meta data columns
snp.meta(kw.dat)    <- snp.meta(kw.dat)[,c(1,2)] # rename snp metadata 
kw.dat <- filter_snps(kw.dat[pop=-"NAPxMON"], maf = 0.05, maf_facets= "pop", remove_garbage = 0.25, min_ind = 0.75, min_loci = 0.75, hf_hets = 0.7) # filter and remove the maternal cross in this dataset 

saveRDS(kw.dat, "kw_rnaseq_standarized.RDS")

## Kellet's whelk GT-Seq unpublished data ## 
kw.gtseq.dat <- readRDS("original_data/kw_gtseq_adults.RDS")

sample.meta(kw.gtseq.dat) <- sample.meta(kw.gtseq.dat)[,c(1,12)]
colnames(sample.meta(kw.gtseq.dat)) <- c("sampID", "pop", ".sample.id")

summarize_facets(kw.gtseq.dat, "pop") # check to see if each pop has more than 5 indivuals 
kw.gtseq.dat <- kw.gtseq.dat[pop=-"CAY"]

kw.gtseq.dat <- filter_snps(kw.gtseq.dat, maf = 0.05, maf_facets= "pop", remove_garbage = 0.25, min_ind = 0.75, min_loci = 0.75, hf_hets = 0.7)

length(unique(sample.meta(kw.gtseq.dat)$pop)) # how many populations 
saveRDS(kw.gtseq.dat, "kw_gtseq_standardized.RDS")

## Monarch data from Hemstrom et al ## 
monarch.dat <- readRDS("original_data/monarch_nomaf.RDS")
sample.meta(monarch.dat) <- sample.meta(monarch.dat)[,1:2]
colnames(sample.meta(monarch.dat)) <- c("sampID", "pop", ".sample.id")
sample.meta(monarch.dat)$sampID <- gsub("\\.", "_", sample.meta(monarch.dat)$sampID)

summarize_facets(monarch.dat, "pop") # check to see if each pop has more than 5 indivuals 
monarch.dat <- monarch.dat[pop = -c("SAI", "VIC")]

monarch.dat <- filter_snps(monarch.dat, maf = 0.05, maf_facets= "pop", remove_garbage = 0.25, min_ind = 0.75, min_loci = 0.75, hf_hets = 0.7) 
length(unique(sample.meta(monarch.dat)$pop)) # how many populations 
saveRDS(monarch.dat, "monarch_standardized.RDS")

## example ranger code for microsatellites
wolf.matrix <- data.table::fread("Yellowstone Montana wolf data set.txt", sep=" ", header=FALSE, colClasses = "character")
wolf.meta <- data.frame(pop=c(rep("Yellowstone", 31), rep("Montana", 66)))
wolf.dat <- snpR:::read_non_biallelic(t(wolf.matrix[,-c(1:2)]), sample.meta = wolf.meta, mDat = "000000") 
sample.meta(wolf.dat)

saveRDS(wolf.dat, "wolf_usat.RDS")


### wolf usat ranger loo
dat <- readRDS("wolf_usat.RDS")
sn <- format_snps(dat, output="pa")                    
train   <- sn[,-c(1)]
sites   <- as.factor(sample.meta(dat)$pop)

rf <- ranger(x=train, importance="permutation", keep.inbag = TRUE, y=sites, num.trees = 10000)
eval <- predict(rf, train, num.trees = 10000)

table(obs=sites, pred=eval$predictions)