library(ranger); library(snpR)
args <- commandArgs(TRUE) 
data <- readRDS(as.character(args[1]))
out <- as.character(args[2])
i <- as.numeric(args[3])

num.trees <- 1000000
set.seed(4+i)

data_sn <- format_snps(data, output="sn")
data_sn <- t(data_sn[,-c(1:2)])
data_sites <- as.factor(sample.meta(data)$pop)

rf_data <- ranger(x=data_sn[-i,], importance="permutation", keep.inbag = TRUE, y=data_sites[-i], num.trees = num.trees, mtry = ncol(data_sn))

# eval <- predict(rf_data, data_sn[i,,drop=FALSE], num.trees = num.trees)

err  <- forestError::quantForestError(rf_data, # the forest 
                                      X.train = data_sn[-i,], # the training data
                                      X.test =  data_sn[i,,drop=FALSE], # the test data
                                      Y.train = data_sites[-i])  # the classifications of the training data

fwrite(data.frame(sample.meta(data)$pop[i], err), file = paste0(out, "_", i, ".txt"), row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")
