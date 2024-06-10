library(snpR)

all_rds <- readRDS("ascertainment_bias/data/monarchs_random.RDS")
setwd("ascertainment_bias/results/structure/all_sites_structure_res/")
clumpp_all <- plot_structure("structure_outfile", k = 4, facet = sample.meta(all_rds)$pop, reps = 10)
file.rename("pop_K4/pop_K4-combined-merged.txt", "pop_K4-combined-merged.txt")

fst_rds <- readRDS("../../../data/monarchs_random_highfst.RDS")
setwd("../high_fst_structure_res/")
clumpp_fst <- plot_structure("structure_outfile", k = 4, facet = sample.meta(all_rds)$pop, reps = 10)

file.rename("pop_K4/pop_K4-combined-merged.txt", "pop_K4-combined-merged.txt")

# to reproduce the plots from the clumpp results, do something like this in the dir with the merged files:
# plot <- plot_structure("merged", k = 4, facet = sample.meta(all_rds)$pop)
