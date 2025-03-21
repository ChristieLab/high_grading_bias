library(snpR)

args <- commandArgs(TRUE)
file <- as.character(args[1])
iter <- as.numeric(args[2])

seed <- 367738 + iter


d <- read_ms(file, chr.length = 10000000)

set.seed(seed)
meta <- data.frame(pop = sample(c("A", "B", "C", "D"), ncol(d), TRUE))
sample.meta(d) <- meta

saveRDS(d, paste0(file, "_converted.RDS"))
