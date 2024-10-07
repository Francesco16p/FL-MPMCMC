## MTM algorithm with random walk proposals
# nsim : length of the chain
# K : number proposed points MTM algorithm 
# log_den : target distribution
# pars : parameters target distribution
# prop_method : covariance structure of the proposal
#               (independent = "ind", antithetic = "anti")
# weight_fun : weighting function MTM 
# sd_0 : initial sd of the proposal distribution 
# alpha_s : optimal acceptance probability (for adaptation)
# d : dimension target density

MTM <- function(nsim, K, log_den ,pars,  prop_method =c("ind","anti"),
                weight_fun = c("identity","sqrt"), sd_0, alpha_s,d)
{
  # Define functions used by the algorithm
  
  require(Matrix)
  
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
  
  ### Weight function MTM ####
  # theta_old : previous value
  # theta_new : proposed point
  # pars : parameters target density
  # weight_fun : weighting function MTM algorithm
  w <- function(theta_new,theta_old, pars, weight_fun, log.o = T) {
    if( weight_fun == "sqrt")
    {
      out <- sqrt( exp( log_den(theta_new, pars)  - log_den(theta_old, pars ))) 
      if(log.o == F)
      {
        return(out)
      }
      return(log(out))
    }
    if( weight_fun == "identity")
    {
      out <- (log_den(theta_new, pars) - log_den(theta_old, pars ) )
      if(log.o == F)
      {
        return(exp(out))
      }
      return(out)
    }
    
  }
  
  ### Proposal distribution for MTM #####################
  # theta_old : current state of the chain
  # tuning_params : tuning parameters for the proposals (see above)
  # prop_method : covariance structure of the proposal
  #               (independent = "ind", antithetic = "anti")
  # nprop_points : number of proposed states MTM algorithm
  # Sigma_cor : covariance proposal distribution (to be implemented in future versions) 
  proposal_mtm <- function(theta_old, tuning_params,prop_method, nprop_points = 1, Sigma_cor =NULL)
  {
    if(prop_method == "ind") 
    {
      d <- length(theta_old)
      sd_vector <- sqrt(tuning_params$lambda*tuning_params$Sigma)
      out <- matrix(rnorm( nprop_points*d, mean = theta_old, sd = sd_vector),
                    ncol = nprop_points, byrow = F )
      return(out)
    }
    if(prop_method  == "anti")
    {
      d <- length(theta_old)
      var_vector <- tuning_params$lambda*tuning_params$Sigma
      if(is.null(Sigma_cor))
      {
        Sigma <- diag(var_vector, ncol = d, nrow = d)
        Psi <- diag(-var_vector/(nprop_points-1+1e-4), ncol = d, nrow = d)
        Sigma_cor <- matrix(0, nrow =  d*nprop_points , ncol = d*nprop_points)
        # Fill Sigma_cor matrix
        for(i in 0:(nprop_points-1))
        {
          for(j in 0:(nprop_points-1))
          {
            if(i == j)  Sigma_cor[d*i+(1:d),d*j+(1:d)] <- Sigma
            else Sigma_cor[d*i+(1:d),d*j+(1:d)] <- Psi
          }
        }
        # Cholesky decomposition
        Sigma_cor_sparse <- as(Sigma_cor, "sparseMatrix")
        U <- chol(Sigma_cor_sparse)
        # Simulation proposed points
        out <-  rep(theta_old,nprop_points) +  t(U)%*%rnorm(d*nprop_points) 
        out <- matrix(out,ncol = nprop_points,byrow = F)
        return(out)
      }
      
    }
  }
  
  # Proposal MTM from theta_new. Needed for computing the final acceptance probability.  
  # theta_old : current state of the chain
  # tuning_params : tuning parameters for the proposals (see above)
  # prop_method : covariance structure of the proposal
  #               (independent = "ind", antithetic = "anti")
  # nprop_points : number of proposed states MTM algorithm
  # Sigma_cor : covariance proposal distribution (to be implemented in future versions) 
  proposal_mtm_2nd_step <- function(theta_new, theta_old, tuning_params,prop_method, nprop_points = 1, Sigma_cor =NULL)
  {
    if(prop_method == "ind") 
    {
      d <- length(theta_new)
      sd_vector <- sqrt(tuning_params$lambda*tuning_params$Sigma)
      out <- matrix(rnorm( (nprop_points)*d, mean = theta_new, sd = sd_vector),
                    ncol = nprop_points, byrow = F )
      return(cbind(theta_old,out))
    }
    if(prop_method  == "anti")
    {
      d <- length(theta_new)
      var_vector <- tuning_params$lambda*tuning_params$Sigma
      if(is.null(Sigma_cor))
      {
        Sigma <- diag(var_vector, ncol = d, nrow = d)
        Psi <- diag(-var_vector/((nprop_points+1)-1+1e-4), ncol = d, nrow = d) # +1 because also theta old is a point
        Sigma_cor <- matrix(0, nrow =  d*(nprop_points+1) , ncol = d*(nprop_points+1))
        # Fill Sigma_cor matrix
        for(i in 0:(nprop_points)) #here we don't have -1 because we are generating the cov assuming also theta_old as point 
        {
          for(j in 0:(nprop_points))
          {
            if(i == j)  Sigma_cor[d*i+(1:d),d*j+(1:d)] <- Sigma
            else Sigma_cor[d*i+(1:d),d*j+(1:d)] <- Psi
          }
        }
        # Here we need to consider that theta_new and theta_old are given
        Sigma_cond <- Sigma_cor[-(1:d),-(1:d)] - Sigma_cor[-(1:d),(1:d)]%*%diag(1/var_vector)%*% Sigma_cor[(1:d),-(1:d)]   # conditional variance
        mu_cond <- rep(theta_new,nprop_points) +  Sigma_cor[-(1:d),(1:d)]%*%diag(1/var_vector)%*%(theta_old - theta_new) # conditional mean
        # Cholesky decomposition
        Sigma_cond_sparse <- as(Sigma_cond, "sparseMatrix")
        U <- chol(Sigma_cond_sparse)
        # Simulation proposed points
        out <-  mu_cond +  t(U)%*%rnorm(d*nprop_points)
        out <- cbind(theta_old, matrix(out,ncol = nprop_points,byrow = F))
        return(out)
      }
      
    }
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
                   pars, prop_method, weight_fun, K)
  {
    
    # Simulations K-1 points from theta_new
    sims_theta_new <- proposal_mtm_2nd_step(theta_new,theta_old, tuning_params,prop_method, nprop_points = K-1)
    
    # Acceptance probability mtm
    log_num_pr <- ( log_den(theta_new,pars) + 
                      w(theta_new= theta_old, theta_old = theta_new,pars, weight_fun,log.o = T)-
                      log( sum( w(theta_new=sims_theta_new,
                                  theta_old = theta_new, pars, weight_fun,log.o = F)
                      ) ) )
    
    log_den_pr <- ( log_den(theta_old,pars) + 
                      w(theta_new= theta_new, theta_old = theta_old, pars, weight_fun,log.o = T)-
                      log( sum( w(theta_new= theta_props,
                                  theta_old = theta_old, pars, weight_fun,log.o = F)
                      ) ) )
    
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
  
  for (i in 2:nsim) {
    # Previous step
    theta_old <- drop(theta_sim[i - 1, ]) # each state of the chain is a row vector
    
    # Propose K values
    theta_props <- proposal_mtm(theta_old, tuning_params, prop_method, nprop_points = K)
    
    # Select one of these proposals with probabilities proportional to the weighting function
    probs1 <-  w(theta_props,theta_old, pars, weight_fun, log.o = F)
    probs1 <- probs1/sum(probs1)
    # Select one of the proposals
    id.sel <- sample(1:K, size = 1, prob = probs1)
    theta_new <- drop(theta_props[,id.sel])
    
    # Accept/Reject the proposed state
    prob_acc[i] <- AMTP(theta_old, theta_new, theta_props, tuning_params,
                        pars, prop_method, weight_fun, K)
    accept[i] <- rbinom(1, 1, prob_acc[i])
    
    if (accept[i] == 1) {
      theta_sim[i, ] <- theta_new
    } else {
      theta_sim[i, ] <- theta_old 
    }
    
    # Adjust probability of acceptance
    tuning_params <- ada_step(tuning_params, drop(theta_sim[i, ]), i^(-0.6), prob_acc[i], alpha_s)
  }
  
  return(list(theta_sim = theta_sim, accept = accept, prob_acc = prob_acc, tuning_params = tuning_params))
}
