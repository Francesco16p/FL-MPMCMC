# First derivative logistic regression 
# theta : value of the parameter on which to evaluate the derivative
# pars : parameters posterior distribution
#        y : response vector
#        z : design matrix
#        hyper_param : sd Gaussian prior 
dlog_logistic <- function(theta, pars)
{
  theta <- as.matrix(theta)
  pr <- 1/( 1 + exp(-pars$z%*%theta) )
  out <- t(pars$z)%*%( pars$y - pr ) - apply(theta,2, function(s) s/(pars$hyper_param^2 ))
  return(out)
}

# Derivative standard multivariate normal
# theta : value of the parameter on which to evaluate the derivative
# pars : not used, retained for compatibility
dlog_stdnorm <- function(theta,pars =NULL)
{
  return(-theta)
}