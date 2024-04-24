monarchs <- readRDS("../assignment_tests/data/monarch_nomaf.RDS")
monarchs <- monarchs[pop = "NAM"]
monarchs <- filter_snps(monarchs, maf = 0.05)
sample.meta(monarchs)$pop <- sample(c("A", "B", "C", "D"), nsamps(monarchs), TRUE)

monarchs <- calc_global_fst(monarchs, "pop")
fst <- get.snpR.stats(monarchs, "pop", "fst")
high_fst <- which(fst$pairwise$fst >= quantile(fst$pairwise$fst, .95, na.rm = TRUE))
high_fst <- snpR:::.paste.by.facet(fst$pairwise, c("group", "position"), "_")[high_fst]
high_fst <- match(high_fst,
                  snpR:::.paste.by.facet(snp.meta(monarchs), c("group", "position"), "_"))


monarchs_high <- monarchs[high_fst,]
monarchs_high <- calc_global_fst(monarchs_high, "pop")
plot_clusters(monarchs_high, "pop")

struc <- plot_structure(monarchs_high, "pop", 
                        method = "structure", 
                        structure_path = "/usr/bin/structure.exe",
                        k = 2:4, iterations = 500, burnin = 100)

plot_clusters(monarchs_high, "pop", "pca")

plot_clusters(monarchs_high, "pop", "tsne")

struc <- plot_structure(monarchs_high, "pop", 
                        method = "snmf", 
                        structure_path = "/usr/bin/structure.exe",
                        k = 2:4, iterations = 500, burnin = 100)

struc <- plot_structure(monarchs_high, "pop", 
                        method = "snapclust", 
                        structure_path = "/usr/bin/structure.exe",
                        k = 2:4, iterations = 500, burnin = 100)
