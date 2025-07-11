library(scrm)
iter <- as.numeric(commandArgs(TRUE)[1])

seed <- 37636 + iter

chrsize <- 10000000
mu <- 1e-8
Ne <- 10000
r <- 2*(1e-8)

gc <- 1000
chrnum <- 10

theta <- 4*Ne*mu*chrsize

rho <- 4*Ne*r*(chrsize - 1)

options(scipen = 999)


call <- paste0(gc, " ", chrnum, " -t ", theta, " -r ", rho, " ", chrsize)
set.seed(seed)
res <- scrm(call, paste0("ms_out_", iter, ".ms"))

