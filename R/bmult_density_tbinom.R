# Evaluates Log Density Of Unnormalized And Unrestricted Beta Distribution
# 
# This function evaluates the log density of the unnormalized and unrestricted Beta distribution.
# It takes probit-transformed samples from a truncated Beta distribution, transforms them back into
# Beta parameter values and evaluates the density
#
#Parameters
#xi_mat matrix with samples from truncated Beta density. These samples should be transformed, so they range
#               over the entire real line
# boundaries list containing indices for upper and lower truncation boundaries
# mult_equal multiplicative elements for each lower and upper bound of each inequality constrained parameter.
# nr_mult_equal vector of multiplicative elements of collapsed parameters
# nr_mult_free vector of multiplicative elements of free parameters
# hyp_direction specifies whether the imposed inequality constrained imposes an increasing (i.e., 'smaller') or
#                      decreasing (i.e., 'larger') trend
#  prior_and_data a numeric vector containing the concentration parameters (when evaluating the prior distribution)
#                       or the updated concentration parameters (when evaluating the posterior distribution)
log_unnormalized_tbinom         <- function(xi_mat, boundaries, binom_equal, hyp_direction, a, b){
  nsamples <- nrow(xi_mat)
  nparam   <- length(a)
  # 1. Transform xi's back to theta's
  backtrans_result   <- tbinom_backtrans(xi_mat, boundaries, binom_equal, hyp_direction)
  theta_mat          <- backtrans_result[['theta_mat']]
  lower_mat          <- backtrans_result[['lower_mat']]
  upper_mat          <- backtrans_result[['upper_mat']]
  
  # 2. Evaluate posterior density for theta's
  logJ_vec  <- rowSums(dnorm(xi_mat, log = TRUE) + log(upper_mat - lower_mat))
  
  if(nparam > 1){
    
    logBeta_vec <- rowSums(t(apply(theta_mat, 1, function(x) dbeta(x, a, b, log = TRUE))))    
    
  } 
  
  if(nparam == 1){
    
    theta_mat   <- as.numeric(theta_mat)
    logBeta_vec <- dbeta(theta_mat, a, b, log = TRUE)  
    
  }
  
  
  logposterior_vec <- logJ_vec + logBeta_vec
  
  return(logposterior_vec)
}