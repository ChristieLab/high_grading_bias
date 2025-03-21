library(coala); library(data.table); library(snpR)

N0 <- 10000 # starting population size (at present)
mu <- 1e-8 # per-base mutation rate
len <- 10000000 # chromosome length
r <- 1 # average number of recombination events per chr
g <- 1 # years per gen
t1 <- 500 # time of ingroup split
m <- 1/N0 # number of individuals moving into each other population per gen
samples <- 50 # sample size
nloci <- 10 # number of chrs
npops <- 8 # number of populations

output_file_name <- "../data/sim_8_islands" # no file extension


# define model
model_basic <- coal_model(sample_size = rep(samples, npops), loci_number = nloci, ploidy = 2, loci_length = len) +
  feat_mutation(rate = par_expr(4*N0*mu*len)) + # 4*N*mu, where mu is mut. rate per locus
  feat_recombination(rate = par_expr(4*N0*r)) # 4*N*r, where r is the probability that a recombination event within the locus occurs in one generation.

for(i in 2:npops){
  model_basic <- model_basic + 
    feat_pop_merge(par_expr(t1/(4*N0*g)), i, 1)
}

for(i in 1:(npops - 1)){
  cat(i, "\n")
  for(j in (i + 1):npops){
    cat("\t", j, "\n")
    model_basic <- model_basic + 
      feat_migration(rate = par_expr(4 * N0 * m), pop_from = i, pop_to = j) +
      feat_migration(rate = par_expr(4 * N0 * m), pop_from = j, pop_to = i)
  }
  cat("=======================\n")
}

model_basic <- model_basic + 
  sumstat_seg_sites() + # generate segregating snps
  par_named("N0") +
  par_named("mu") +
  par_named("r") +
  par_named("g") +
  par_named("t1") +
  par_named("m")

set.seed(123)
sim <- simulate(model_basic, seed = 123, pars = c(N0 = N0, len = len,
                                                  r = r, g = g, t1 = t1, 
                                                  mu = mu, m = m))

saveRDS(sim, paste0(output_file_name, "raw.RDS"))


# parse
genos <- sim$seg_sites
meta <- vector("list", length(genos))
for(i in 1:length(genos)){
  meta[[i]] <- data.table(chr = i, position = floor(as.numeric(genos[[i]]$position)*len))
  genos[[i]] <- genos[[i]]$snps[seq(1, nrow(genos[[i]]), 2),] +
    genos[[i]]$snps[seq(2, nrow(genos[[i]]), 2),]
  genos[[i]] <- t(genos[[i]])
  genos[[i]] <- as.data.table(genos[[i]])
}

meta <- rbindlist(meta)
genos <- rbindlist(genos)


# send to vcf
srd <- import.snpR.data(genos, sample.meta = data.frame(pop = rep(LETTERS[1:npops], each = 50)),
                        snp.meta = meta)


saveRDS(srd, paste0(output_file_name, "snpR.RDS"))
