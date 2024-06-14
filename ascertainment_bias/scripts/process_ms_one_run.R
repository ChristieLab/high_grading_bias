args <- commandArgs(TRUE)
dat <- as.character(args[1])
iter <- as.numeric(args[2])

n_boots <- 100

i <- 1

cmd <- "squeue -u whemstro -h -t pending,running -r | wc -l"

cat("Starting queuing process for: ", dat, ", which is iter:", iter, "\n")

# queue up each bootstrap set for this iter, checking to see if too many jobs are queued already before queueing
wait_status <- 1
while(i <= n_boots){
  num_running <- as.numeric(system(cmd, intern = TRUE))
  
  cat("Attempt: ", wait_status, ", ")
  if(num_running < 990){
    cat(num_running, "jobs running, queuing boot ", i, ".\n")
    tcmd <- paste0("sbatch process_ms_one_run.sh ",
                   dat, " ", iter, " ", i)
    system(tcmd)
    i <- i + 1
    wait_status <- 1
  }
  else{
    cat(num_running, "jobs running, waiting...\n")
    wait_status <- wait_status + 1
    Sys.sleep(10) # check back in a while
  }
}
