library(GeneArchEst); library(data.table)
args <- commandArgs(TRUE)
run <- as.numeric(args[1])
outfile <- as.character(args[2])

d <- process_ms("../data/theta4k_1000_10_rho40k.txt", 10000000)
genotypes <- list(d$x[,1:250], d$x[, 251:500],
                  d$x[,501:750], d$x[,751:1000])

set.seed(123344 + run)
effects <- rep(0, nrow(d$meta))
meta <- d$meta

effects <- cbind(effects, effects, effects, effects)

rate <- .05/3 # pev: none
migration <- matrix(rate, 4, 4) # migration matrix, from row into column
diag(migration) <- 1-(rate*3)
migration

# sim for 50 gens
set.seed(123345)
tno <- gs(genotypes, meta = meta, effects = effects, h = .5, gens = 50, chr.length = 10000000,
          migration = migration,
          growth.function = function(n) BL_growth(n, 2), # use negative binomial instead
          mutation = 1e-8,
          var.theta = c(0.5, 0.5, 0.5, 0.5),
          mutation.effect.function = function(n) rep(0, n),
          K_thin_post_surv = 125,
          thin = FALSE,
          thin_fixed = TRUE,
          survival.function = function(phenotypes, opt_pheno) rep(1, length(phenotypes)),
          selection.shift.function = function(opt, ...) opt,
          plot_during_progress = TRUE,
          verbose = TRUE)
saveRDS(tno, paste0(outfile, "_raw.RDS"))

cat("=================================================\nWARNINGS:\n\n")
print(warnings())
cat("\n\nEND WARNINGS\n=================================================\n")

library(snpR)
genos <- cbind(tno$genotypes[[1]], tno$genotypes[[2]], tno$genotypes[[3]], tno$genotypes[[4]])
genos <- convert_2_to_1_column(genos)
genos <- t(genos)
d <- import.snpR.data(as.data.frame(genos), snp.meta = tno$meta, 
                      sample.meta = data.frame(pop = c(rep("A", ncol(tno$genotypes[[1]])/2),
                                                       rep("B", ncol(tno$genotypes[[2]])/2),
                                                       rep("C", ncol(tno$genotypes[[3]])/2),
                                                       rep("D", ncol(tno$genotypes[[4]])/2))))
saveRDS(d, paste0(outfile, ".RDS"))
