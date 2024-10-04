### Log-densities ###

### Posterior logistic regression ###
# theta : value of the parameter on which to evaluate the derivative
# pars : list containing the parameters posterior distribution 
#        y : response vector
#        z : design matrix
#        hyper_param : sd Gaussian prior 
log_post_logistic <- function(theta, pars) {
  theta <- as.matrix(theta)
  y <- pars$y
  z <- pars$z
  sd_p <- pars$hyper_param
  # parameter dimension
  d <- dim(z)[2]
  # linear predictors
  eta <- as.matrix(z%*%theta)
  # log-likelihood and log prior
  log_lik <- apply(eta, 2, function(s)  sum(y*s - log1p( exp(s) )  ) ) 
  log_prior <- apply(theta, 2, function(s) sum( -0.5*(s/sd_p)^2 ) )
  log_post <- log_prior + log_lik
  # The output is a row vector
  return(log_post)
}

# Un-normalized multivariate Gaussian
# theta : value of the parameter on which to evaluate the derivative
# pars : list containing the parameters of the distribution
log_Gauss <- function(theta, pars)
{
  mu <- pars$mu
  prec <- pars$prec
  if(is.matrix(theta))
  {
    if( dim(theta)[2]>1)
    { 
      out <- rep(NA,dim(theta)[2])
      for(j in 1:dim(theta)[2])
      {
        out[j] <- -0.5*sum( theta[,j]*(prec%*%theta[,j]))
      }
    }  
    else{
        out <- -0.5*sum( theta*(prec%*%theta) )
      }
  }
  else{
    out <- -0.5*sum( theta*(prec%*%theta) )
  }
  
  return(out)
}

### Un-normalized multivariate Gaussian ###
# theta : value of the parameter on which to evaluate the derivative
# pars : not used, retained for compatibility
log_Gauss <- function(theta, pars)
{
  mu <- pars$mu
  prec <- pars$prec
  if(is.matrix(theta))
  {
    if( dim(theta)[2]>1)
    { 
      out <- rep(NA,dim(theta)[2])
      for(j in 1:dim(theta)[2])
      {
        out[j] <- -0.5*sum( theta[,j]*(prec%*%theta[,j]))
      }
    }  
    else{
      out <- -0.5*sum( theta*(prec%*%theta) )
    }
  }
  else{
    out <- -0.5*sum( theta*(prec%*%theta) )
  }
  
  return(out)
}
