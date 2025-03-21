### RUBIAS loo self assignment ====
library(rubias); library(tidyverse); library(snpR)
        
# setwd("~/assignment_power_analyses/data")
args <- commandArgs(TRUE) 
data <- readRDS(as.character(args[1]))
outfile <- as.character(args[2])

### format input data 

### Create functions to sort data and run the sa function in rubias 
### the data frame must have these 4 columns: sample_type, repunit, collection, indiv
rubias_data_setup <- function(data){
  sample.meta(data)$indiv              <- sample.meta(data)[,1]
  # collection populations
  data_rafm <- format_snps(data, facets="indiv", "rafm") 
  data_rafm$collection         <- as.character(sample.meta(data)[,2])
  data_rafm$sample_type        <- "reference"
  data_rafm$repunit            <- "NA" # NA as we are just using the
  data_rafm <- data_rafm %>%
    select(indiv, sample_type, repunit, collection, everything())
  return(data_rafm)
  }
run_rubias_sa <- function(data_rafm){
  sa <- self_assign(reference = data_rafm, gen_start_col = 5)
  sa_to_site <- sa %>%
    group_by(indiv, collection, repunit, inferred_collection) %>%
    summarise(repu_scaled_like = sum(scaled_likelihood)) %>% 
    group_by(indiv) %>%
    filter(repu_scaled_like==max(repu_scaled_like)) 
  return(sa_to_site)
  # add line to write result text files and include dataset name
}

### Run the codes 
rubias_input <- rubias_data_setup(data)
sa_result    <- run_rubias_sa(rubias_input)

##   sa_to_site   <- sa_result %>% 
##     #group_by(indiv, collection, repunit, inferred_collection) %>%
##     #summarise(repu_scaled_like = sum(scaled_likelihood)) %>% 
##     group_by(indiv) %>%
##     filter(repu_scaled_like==max(repu_scaled_like)) %>%
##     summarise(percentage_match = sum(collection == inferred_collection) / n()) 
##   mean(sa_to_site$percentage_match)

#=========save=======================
outfile <- paste0(outfile, ".txt")
fwrite(sa_result, outfile, col.names = TRUE, row.names = FALSE, sep = "\t")
