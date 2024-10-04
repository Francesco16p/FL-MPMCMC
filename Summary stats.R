### Summary statistics for simulation results ###

### Expected squared jump distance ###
# sim : vector containing the simulated chain
ESJD <- function(sim) {
  out <- mean((sim[-1] - sim[-length(sim)])^2)
  return(out)
}

### Output simulation ###
# sim : Matrix containing the simulated chain.
#       Different rows represent different iterations
summary_stats <- function(sim)
{
  require(coda)
  t <- dim(sim)[1] 
  sim <- sim[round(t/4):t,] # remove burn-in
  ESS_sum <- effectiveSize(as.mcmc(sim))
  ESJD_sum <- apply(sim,2,ESJD)
  IAT_sum <- dim(sim)/ESS_sum 
  return(list(ESS_sum = ESS_sum,ESJD_sum = ESJD_sum, IAT_sum =  IAT_sum))
}
