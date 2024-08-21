library(GeneArchEst); library(data.table)

args <- commandArgs(TRUE)
parms <- as.matrix(data.table::fread(as.character(args[1])))
iter <- as.numeric(args[2])
outfile <- as.character(args[3])
num_effect <- unlist(as.numeric(parms[iter,1]))
mig_rate <- unlist(as.numeric(parms[iter,2]))
optima <- unlist(as.numeric(parms[iter,3:6]))
omega <- unlist(as.numeric(parms[iter,7]))
num_offspring <- unlist(as.numeric(parms[iter,8]))
K <- unlist(as.numeric(parms[iter,9]))
seed <- unlist(as.numeric(parms[iter, 10]))



d <- process_ms("../data/theta4k_1000_10_rho40k.txt", 10000000)
genotypes <- list(d$x[,1:250], d$x[, 251:500],
                  d$x[,501:750], d$x[,751:1000])

set.seed(123343 + seed)
prop_effect <- num_effect/nrow(d$x)

effects <- rbayesB_fixed(nrow(d$x), num_effect, 5, 1)
meta <- d$meta

effects <- cbind(effects, effects, effects, effects)

rate <- mig_rate/3
migration <- matrix(rate, 4, 4) # migration matrix, from row into column
diag(migration) <- 1-(rate*3)
migration

#==============with selection=====================
# sim for 50 gens
set.seed(123344 + seed)
#tno <- gs(genotypes, meta = meta, effects = effects, h = .5, gens = 50, chr.length = 10000000,
#          migration = migration,
#          growth.function = function(n) BL_growth(n, num_offspring), # use negative binomial instead
#          mutation = 1e-8,
#          var.theta = c(0.5, 0.5, 0.5, 0.5),
#          mutation.effect.function = function(n) rbayesB(n, 1-prop_effect, 5, 1),
#          K_thin_post_surv = K,
#          thin = FALSE,
#          thin_fixed = TRUE,
#          starting.surv.opt = optima,
#          survival.function = function(phenotypes, opt_pheno) BL_survival(phenotypes, opt_pheno, omega),
#          selection.shift.function = function(opt, ...) opt,
#          plot_during_progress = FALSE, 
#          sampling_point = "parents",
#          verbose = TRUE)

#saveRDS(tno, paste0(outfile, "_raw.RDS"))

#if("genotypes" %in% names(tno)){ 
#  library(snpR)
#  genos <- cbind(tno$genotypes[[1]], tno$genotypes[[2]], tno$genotypes[[3]], tno$genotypes[[4]])
#  genos <- convert_2_to_1_column(genos)
#  genos <- t(genos)
#  d <- import.snpR.data(as.data.frame(genos), snp.meta = tno$meta[[1]],
#                        sample.meta = data.frame(pop = c(rep("A", ncol(tno$genotypes[[1]])/2),
#                                                         rep("B", ncol(tno$genotypes[[2]])/2),
#                                                         rep("C", ncol(tno$genotypes[[3]])/2),
#                                                         rep("D", ncol(tno$genotypes[[4]])/2))))
#  saveRDS(d, paste0(outfile, ".RDS"))
#}

#rm(tno)

#==============control=====================
# sim for 50 gens
effects <- rep(0, nrow(d$x))
effects <- cbind(effects, effects, effects, effects)

set.seed(123344 + seed)
tno <- gs(genotypes, meta = meta, effects = effects, h = .5, gens = 50, chr.length = 10000000,
          migration = migration,
          growth.function = function(n) BL_growth(n, num_offspring), # use negative binomial instead
          mutation = 1e-8,
          var.theta = c(0.5, 0.5, 0.5, 0.5),
          mutation.effect.function = function(n) rep(0, n),
          K_thin_post_surv = K,
          thin = FALSE,
          thin_fixed = TRUE,
          starting.surv.opt = optima,
          survival.function = function(phenotypes, opt_pheno) BL_survival(phenotypes, opt_pheno, omega),
          selection.shift.function = function(opt, ...) opt,
          plot_during_progress = FALSE, 
          sampling_point = "parents",
          verbose = TRUE)

#saveRDS(tno, paste0(outfile, "_raw.RDS"))

if("genotypes" %in% names(tno)){ 
  library(snpR)
  genos <- cbind(tno$genotypes[[1]], tno$genotypes[[2]], tno$genotypes[[3]], tno$genotypes[[4]])
  genos <- convert_2_to_1_column(genos)
  genos <- t(genos)
  d <- import.snpR.data(as.data.frame(genos), snp.meta = tno$meta[[1]],
                        sample.meta = data.frame(pop = c(rep("A", ncol(tno$genotypes[[1]])/2),
                                                         rep("B", ncol(tno$genotypes[[2]])/2),
                                                         rep("C", ncol(tno$genotypes[[3]])/2),
                                                         rep("D", ncol(tno$genotypes[[4]])/2))))
  saveRDS(d, paste0(outfile, "_control.RDS"))
}
