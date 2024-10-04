## MTM algorithm with MALA proposal
# nsim : length of the chain
# K : number proposed points MTM algorithm 
# log_den : target distribution
# dlog_den : first derivative target distribution
# pars : parameters target distribution
# weight_fun : weighting function MTM
# sd_0 : initial sd of the proposal distribution 
# alpha_s : optimal acceptance probability (for adaptation)
# d : dimension target density
# adaptive : (logical) should adaptation of the proposal be used?

MTM_MALA <- function(nsim, K, log_den , dlog_den, pars,
                weight_fun = c("identity","sqrt"), sd_0, alpha_s,d, adaptive = TRUE)
{
  # Define functions used by the algorithm
  
  ### log-sum exp function ###
  logSumExp <- function(values)
  {
    m <- max(values)
    return(m + log( sum( exp(values - m) ) ))
  }
  
  ###  Adaptive step-size algorithm ###
  # tuning_params : parameters at time t
  #                 - lambda : scale parameter
  #                 - Sigma : covariance matrix 
  #                 - mu : mean, only used to estimate Sigma at time t + 1
  # theta_new : current state of the chain
  # m : learning rate
  # alpha : acceptance probability time t+1
  # alpha_s : optimal alpha
  ada_step <- function(tuning_params, theta_new, m, alpha, alpha_s) {
    log_lambda <- log(tuning_params$lambda)
    Sigma <- tuning_params$Sigma
    mu <- tuning_params$mu
    d <- length(mu)
    #update log_lambda
    log_lambda <- log_lambda + m * (alpha - alpha_s)
    diff_theta_mu <- theta_new - mu
    # update vector sd (not used in this version)
    #  Sigma <-  Sigma + m*(diff_theta_mu^2- Sigma)
    #  mu <- mu + m*diff_theta_mu
    return(list( lambda =exp(log_lambda), mu = mu, Sigma = Sigma ))
  }
  
  ### log density MALA proposal ###
  # theta_new : proposed theta
  # theta_old : previous value for theta
  # tuning_params : tuning parameters for the proposals (see above)
  log_mala <- function(theta_new,theta_old,tuning_params)
  {
    dim_old <- dim(as.matrix(theta_old))[2]
    dim_new <- dim(as.matrix(theta_new))[2]
    if(dim_old==1)
    {
      sd_vector <- sqrt(tuning_params$lambda*tuning_params$Sigma)
      m <- theta_old + (sd_vector^2)*dlog_den(theta_old,pars)/2
      out <- apply(as.matrix(theta_new),2, function(s) -0.5*sum(  ((s-m)/sd_vector)^2 ) )
      return(out)
    }
    if( (dim_old > 1)&(dim_new == 1) )  
    {
      sd_vector <- sqrt(tuning_params$lambda*tuning_params$Sigma)
      m <- theta_old + (sd_vector^2)*dlog_den(theta_old,pars)/2
      out <- apply( m ,2, function(s) -0.5*sum(  ((theta_new-s)/sd_vector)^2 ) )
      return(out)
    }
    else
    {
      stop("conditions not verified in function log_mala")
    }
   
  }

  ### Weighting function MTM ####
  # theta_old : previous value
  # theta_new : proposed point
  # pars : parameters target density
  # weight_fun : weighting function MTM algorithm
  # tuning_params : tuning parameters for the proposals (see above) 
  w <- function(theta_new,theta_old, pars, weight_fun,tuning_params) {
    if( weight_fun == "sqrt")
    {
      out <- (0.5) * (log_den(theta_new,pars) - log_den(theta_old,pars)  
                        + log_mala(theta_new = theta_old, theta_old=theta_new, tuning_params)
                        - log_mala(theta_new = theta_new, theta_old=theta_old, tuning_params) ) 
      return(out)
    }
    if( weight_fun == "identity")
    {
      out <- (log_den(theta_new,pars) - log_den(theta_old,pars)  
              + log_mala(theta_new = theta_old, theta_old=theta_new, tuning_params)
              - log_mala(theta_new = theta_new, theta_old=theta_old, tuning_params) )
      return(out)
    }
    
  }
  
  ### Proposal distribution for MTM MALA ###
  # theta_old : current state of the chain
  # tuning_params : tuning parameters for the proposals (see above)
  # nprop_points : number of proposed states MTM algorithm
  proposal_mtm <- function(theta_old, tuning_params, nprop_points)
  {
    sd_vector <- sqrt(tuning_params$lambda*tuning_params$Sigma)
    m <- theta_old + (sd_vector)^2*dlog_den(theta_old,pars)/2
    out <- ( matrix(rep(m,nprop_points), ncol = nprop_points, byrow = F) +
             matrix( rnorm(d*nprop_points, sd = sd_vector), ncol = nprop_points, byrow = F ) )
    return(out)
  }
  
  ### Probability of acceptance for a given set of proposals MTM ###
  # theta_old : previous value of the chain
  # theta_new : point selected among the K proposed
  # theta_props : proposed states  (must include theta_new)
  # tuning_params : tuning parameters for the proposals (see above)
  # pars : parameters target distribution
  # weight_fun : weighting function MTM
  # K : number of proposed points
  AMTP <- function(theta_old, theta_new, theta_props, tuning_params,
                   pars, weight_fun, K)
  {
    
    # Simulations K-1 points from theta_new
    sims_theta_new <- proposal_mtm(theta_new, tuning_params, K-1)
    sims_theta_new <- cbind(theta_old,sims_theta_new)
    
    # Acceptance probability mtm
    log_num_pr <-( log_den(theta_new,pars) + log_mala(theta_new = theta_old, theta_old =theta_new, tuning_params )+
                   w(theta_new = theta_old,theta_old = theta_new,pars = pars,weight_fun = weight_fun,
                     tuning_params = tuning_params) -
                  logSumExp( 
                    w(theta_new = sims_theta_new,theta_old = theta_new,pars = pars,weight_fun = weight_fun,
                      tuning_params = tuning_params)
                  ))
     
    log_den_pr <- ( log_den(theta_old,pars) + log_mala(theta_new = theta_new, theta_old =theta_old, tuning_params )+
                      w(theta_new = theta_new,theta_old = theta_old,pars = pars,weight_fun = weight_fun,
                        tuning_params = tuning_params) -
                      logSumExp( 
                        w(theta_new = theta_props,theta_old = theta_old, pars = pars,weight_fun = weight_fun,
                          tuning_params = tuning_params)
                      ))
    
    pr_acc <- min( 1, exp(log_num_pr - log_den_pr))
    if(is.na(pr_acc))
    {
      warning("NA acceptance probability in MTM")
      return(0)
    }
    return(pr_acc)
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
  log_lambda <- rep(NA,nsim)
  #theta_sim[1, ] <- runif(d,-2,2)
  
  for (i in 2:nsim) {
    # Previous step
    theta_old <- drop(theta_sim[i - 1, ]) # each state of the chain is a row vector
    
    # Propose K values
    theta_props <- proposal_mtm(theta_old, tuning_params, nprop_points = K)
    # Select one of these proposals with probabilities proportional to the weighting function
    probs1 <-  w(theta_new =  theta_props,theta_old = theta_old, pars, weight_fun,tuning_params)
    probs1 <- exp(probs1 - logSumExp(probs1))
    # Select one of the proposals
    if(any(is.na(probs1))) browser()
    id.sel <- sample(1:K, size = 1, prob = probs1)
    theta_new <- drop(theta_props[,id.sel])
    # Accept/Reject the proposed state
    prob_acc[i] <- AMTP(theta_old, theta_new, theta_props, tuning_params,
                        pars, weight_fun, K)
    accept[i] <- rbinom(1, 1, prob_acc[i])
    
    if (accept[i] == 1) {
      theta_sim[i, ] <- theta_new
    } else {
      theta_sim[i, ] <- theta_old 
    }
    
    # Adjust probability of acceptance
    if(adaptive == T)
    {
      tuning_params <- ada_step(tuning_params, drop(theta_sim[i, ]), i^(-0.6), prob_acc[i], alpha_s)
    }
    log_lambda[i] <- log(tuning_params$lambda)
  }
  return(list(theta_sim = theta_sim, accept = accept, prob_acc = prob_acc, tuning_params = tuning_params,
              log_lambda = log_lambda))
}