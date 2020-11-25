# Evaluates Log Density Of Unnormalized And Unrestricted Dirichlet Distribution
# 
# This function evaluates the log density of the unnormalized and unrestricted Dirichlet distribution.
# It takes probit-transformed samples from a truncated Dirichlet distribution, transforms them back into
# Dirichlet parameter values and evaluates the density
#
# @param xi_mat matrix with samples from truncated Dirichlet density. These samples should be transformed, so they range
#               over the entire real line
# @param boundaries list containing indices for upper and lower truncation boundaries
# @param mult_equal multiplicative elements for each lower and upper bound of each inequality constrained parameter.
# @param nr_mult_equal vector of multiplicative elements of collapsed parameters
# @param nr_mult_free vector of multiplicative elements of free parameters
# @param hyp_direction specifies whether the imposed inequality constrained imposes an increasing (i.e., 'smaller') or
#                      decreasing (i.e., 'larger') trend
# @param prior_and_data a numeric vector containing the concentration parameters (when evaluating the prior distribution)
#                       or the updated concentration parameters (when evaluating the posterior distribution)
#
log_unnormalized_tdir         <- function(xi_mat, boundaries, mult_equal, nr_mult_equal, nr_mult_free, hyp_direction, prior_and_data){
  nsamples <- nrow(xi_mat)
  nparam   <- length(prior_and_data)
  # 1. Transform xi's back to theta's
  backtrans_result   <- tdir_backtrans(xi_mat, boundaries, mult_equal, nr_mult_equal, nr_mult_free, hyp_direction)
  theta_mat          <- backtrans_result[['theta_mat']]
  lower_mat          <- backtrans_result[['lower_mat']]
  upper_mat          <- backtrans_result[['upper_mat']]
  
  prior_and_data_mat <- matrix(prior_and_data,
                               ncol  = nparam,
                               nrow  = nsamples,
                               byrow = TRUE)
  
  # 2. Evaluate posterior density for theta's
  logJ_vec  <- rowSums(dnorm(xi_mat, log = TRUE) + log(upper_mat - lower_mat))
  
  if(nparam > 2){
    
    # unnormalized truncated Dirichlet
    logdD_vec <- lgamma(sum(prior_and_data)) - sum(lgamma(prior_and_data)) +
      (rowSums((prior_and_data_mat - 1) * log(theta_mat)))
    
  } else {
    
    # truncated beta
    logdD_vec <- dbeta(theta_mat[,1], prior_and_data[1], prior_and_data[2], log = TRUE)
    
  }
  
  
  logposterior_vec <- logJ_vec + logdD_vec
  
  return(logposterior_vec)
}