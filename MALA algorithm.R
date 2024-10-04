##  Single-proposal MALA algorithm with adaptive proposal
# nsim : length of the chain
# log_den : target distribution
# dlog_den : first derivative target distribution
# pars : parameters target distribution
# sd_0 : initial sd of the proposal distribution 
# alpha_s : optimal acceptance probability (for adaptation)
# d : dimension target density

MALA <- function(nsim,  log_den, dlog_den, pars, sd_0, alpha_s,d)
{
  # Define functions used by the algorithm
  
  ####  Adaptive step-size algorithm #####
  # tuning_params : parameters at time t
  #                 - lambda : scale parameter
  #                 - Sigma : covariance matrix 
  #                 - mu : mean, only used to estimate Sigma at time t + 1
  # m : learning rate
  # alpha : acceptance probability time t+1
  # alpha_s : optimal alpha
  # theta_new : current state of the chain
  
  ada_step <- function(tuning_params, theta_new, m, alpha, alpha_s) {
    log_lambda <- log(tuning_params$lambda)
    Sigma <- tuning_params$Sigma
    mu <- tuning_params$mu
    d <- length(mu)
    #update
    log_lambda <- log_lambda + m * (alpha - alpha_s)
    diff_theta_mu <- theta_new - mu
    # update vector sd
    Sigma <-  Sigma + m*(diff_theta_mu^2- Sigma)
    mu <- mu + m*diff_theta_mu
    
    return(list( lambda =exp(log_lambda), mu = mu, Sigma = Sigma ))
  }
  
  ### Proposal distribution for MTM MALA ###
  # theta_old : current state of the chain
  # tuning_params : tuning parameters for the proposals (see above)
  # pars : parameters target distribution
  proposal_mala <- function(theta_old, tuning_params, pars)
  {
    delta <- tuning_params$lambda*tuning_params$Sigma
    out <- mvtnorm::rmvnorm(1, mean = theta_old + delta*dlog_den(theta_old, pars)/2, sigma = diag(delta) ) 
    return(out_sim = c(out) )
  }
  
  # Allocate the matrix for the simulation
  theta_sim <- matrix(0, nrow = nsim, ncol = d)
  
  # Define starting values for adaptation
  mu <- rep(0,d)
  Sigma <- rep(1,d)
  tuning_params <- list(lambda = sd_0, mu = mu, Sigma = Sigma)
  
  
  # Allocate vectors for tracking the acceptance rate
  accept <- rep(1, nsim)
  prob_acc <- rep(0, nsim)
  
  for (i in 2:nsim) {
    # Previous step
    theta_old <- drop(theta_sim[i - 1, ]) # each state of the chain is a row vector
    
    # Propose a new value
    delta <- tuning_params$lambda*tuning_params$Sigma
    theta_new <- proposal_mala(theta_old, tuning_params, pars)
    
    ### Acceptance probability mtm
   
    # derivative at theta_old
    d_old <- dlog_den(theta_old, pars)
    d_new <- dlog_den(theta_new, pars)
    
    log_num_pr <- ( log_den(theta_new,pars) + mvtnorm::dmvnorm(theta_old, mean = theta_new + d_new*delta/2, 
                                                               sigma = diag(delta),log = TRUE)
                    )
    
    log_den_pr <- ( log_den(theta_old,pars) + mvtnorm::dmvnorm(theta_new, mean = theta_old + d_old*delta/2, 
                                                               sigma = diag(delta),log = TRUE)
                      )
    
    prob_acc[i] <- min( 1, exp(log_num_pr - log_den_pr))
    
    accept[i] <- rbinom(1, 1, prob_acc[i])
    if (accept[i] == 1) {
      theta_sim[i, ] <- theta_new
    } 
    
    else {
      theta_sim[i, ] <- theta_old 
    }
    
    # Adjust probability of acceptance
    tuning_params <- ada_step(tuning_params, drop(theta_sim[i, ]), i^(-0.6), prob_acc[i], alpha_s)
  }
  
  return(list(theta_sim = theta_sim, accept = accept, prob_acc = prob_acc, tuning_params = tuning_params))
}