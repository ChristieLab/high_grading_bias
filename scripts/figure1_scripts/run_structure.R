library(snpR); library(data.table)

args <- commandArgs(TRUE)
parms <- read.table(as.character(args[1]), header = FALSE, sep = "\t")
i <- args[2]
dat <- as.character(args[3])
target_dir <- as.character(args[4])

ns <- 5000

niter <- 100000
burnin <- 20000

dat <- readRDS(dat)

if(!dir.exists(target_dir)){dir.create(target_dir)}

setwd(target_dir)

ndir <- paste0("r", i)
dir.create(ndir)

seed <- parms[i,2]
k <- parms[i,1]

set.seed(seed)
if(nrow(dat) > ns){
  dat <- dat[sample(nrow(dat), ns, replace = FALSE),]
}


setwd(ndir)

plot_structure(dat, method="structure", k=k, 
               structure_path = "/home/hemstrow/bin/structure",
               facet = "pop", 
               iteration = niter, 
               burnin = burnin, 
               reps = 1, 
               clumpp = FALSE, 
               cleanup = FALSE)
