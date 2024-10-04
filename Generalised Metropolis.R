## Generalised Metropolis algorithm 

# nsim = number of iterations required
# K = number of proposed states at each iteration
# log_den = target density
# pars = parameters target density
# sd0 = initial value variance proposal

GMH <- function(nsim, K, log_den ,pars, sd_0, alpha_s,d)
{
  require(Matrix)
  ### log-sum exp function ###
  logSumExp <- function(values)
  {
    m <- max(values)
    return(m + log( sum( exp(values - m) ) ))
  }
  ###  Adaptive step-size algorithm ###
  # theta_new : current state of the chain
  # tuning_params : parameters at time t
  #                 - lambda : scale parameter
  #                 - Sigma : covariance matrix 
  #                 - mu : mean, only used to estimate Sigma at time t + 1
  # m : learning rate
  # alpha : acceptance probability time t+1
  # alpha_s : optimal alpha
  ada_step_GM <- function(tuning_params, theta_new, m, alpha, alpha_s) {
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
  
  ### Proposal distribution for GM ###
  # theta_old : current state of the chain
  # tuning_params : tuning parameters for the proposals (see above)
  # nprop_points : number of proposed states MTM algorithm
  proposal_mtm_GM <- function(theta_old, tuning_params, nprop_points = 1)
  {
    d <- length(theta_old)
    sd_vector <- sqrt(tuning_params$lambda*tuning_params$Sigma)
    rho <- rnorm( d, mean = theta_old, sd = sd_vector/sqrt(2))
    out <- matrix(rnorm( nprop_points*d, mean = rho, sd = sd_vector/sqrt(2)),
                  ncol = nprop_points, byrow = F )
    return(out)
  }
  
  ### Transition matrix generalized metropolis ###
  # theta_props_kp1: matrix with all the K+1 proposed points
  # tuning_params : tuning parameters for the proposals (see above)
  # K : number of proposed states MTM algorithm
  w_gm <- function(theta_props_kp1,tuning_params,K)
  {
    probs <- numeric(K+1)
    
    # log-density of the proposed points
    log_pi <- log_den(theta_props_kp1, pars)
    
    # Transition probabilities
    probs <- exp(log_pi - logSumExp(log_pi))
    return(probs)
  }
  
  # Allocate the matrix for the simulation
  theta_sim <- matrix(0, nrow = nsim, ncol = d)
  
  # Define starting values for adaptation
  mu <- rep(0,d)
  Sigma <- rep(1,d)
  tuning_params <- list(lambda = sd_0, mu = mu, Sigma = Sigma)
  
  
  # Allocate vectors for tracking the acceptance rate
  accept <- rep(0, nsim)
  prob_acc <- rep(0, nsim)
  
  # Indicator variable of the current state
  nu <- 1
  
  for (i in 2:nsim) {
    # Previous step
    theta_old <- drop(theta_sim[i - 1, ]) # each state of the chain is a row vector
    
    # Propose K values
    theta_props <- proposal_mtm_GM(theta_old, tuning_params, nprop_points = K)
    Z_prop <- matrix(0,nrow = d, ncol = K +1)
    Z_prop[,nu] <- theta_old
    Z_prop[,-nu] <- theta_props
    
    # update value of nu
    pr <- w_gm(Z_prop, dataset, hyper_param, tuning_params,K)
    prob_acc[i] <- sum(pr[-nu])
    # Sample a new nu
    nu <- sample(1:(K+1),1,prob = pr) # Simplified acceptance probabilities of Holbrook
    theta_new <- drop(Z_prop[,nu])
    # Accept/Reject the proposed state
    
    accept[i] <- any(theta_new != theta_old)
    theta_sim[i, ] <- theta_new
    
    # Adjust probability of acceptance
    tuning_params <- ada_step_GM(tuning_params, drop(theta_sim[i, ]), i^(-0.6), prob_acc[i], alpha_s)
  }
  
  return(list(theta_sim = theta_sim, accept = accept, prob_acc = prob_acc, tuning_params = tuning_params))
}



