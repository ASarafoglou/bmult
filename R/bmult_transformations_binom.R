#' Transforms Truncated Beta Samples To Real Line
#'
#' Transforms samples from a truncated beta density to the real line using a stick-breaking algorithm.
#' This algorithm is suitable for mixtures of equality constrained parameters, inequality constrained parameters,
#' and free parameters
#'
#' @param theta_mat matrix with samples from truncated beta density
#' @param boundaries list containing indices for upper and lower truncation boundaries
#' @param binom_equal multiplicative elements for each lower and upper bound of each inequality constrained parameter.
#' @param hyp_direction specifies whether the imposed inequality constrained imposes an increasing (i.e., 'smaller') or
#'                      decreasing (i.e., 'larger') trend
#' @return matrix with transformed samples
#' @export
tbinom_trans           <- function(theta_mat, boundaries, binom_equal, hyp_direction){
  # theta_mat : prior or posterior samples (nsamples x nparam)
  # boundaries: list with upper and lower bounds
  nparameters   <- length(boundaries)
  nsamples      <- nrow(theta_mat)
  xi_mat        <- matrix(NA, ncol = nparameters, nrow = nsamples)
  
  # To transform the variables we need to move from the smallest to the largest parameter.
  # Therefore, we need to determine the direction of the hypothesis and adjust the algorithm
  # accordingly.
  if(hyp_direction == 'smaller'){
    
    direction          <- 1:nparameters # move from smallest to largest parameter
    index              <- 1:nparameters # index in xi_mat matrix
    smallest_parameter <- 1
    
  } else {
    
    direction          <- nparameters:1   # move from smallest to largest parameter
    index              <- 1:nparameters   # index in xi_mat matrix
    smallest_parameter <- nparameters
    
  }
  
  for(k in direction){
    
    smaller_values <- boundaries[[k]]$lower
    larger_values  <- boundaries[[k]]$upper
    
    lower     <- rep(0, nsamples)
    
    if(length(smaller_values) != 0){
      
      #binom_mat_lower  <- matrix(binom_equal[[k]]$lower, ncol=length(binom_equal[[k]]$lower), nrow=nrow(theta_mat), byrow=TRUE)
      #lower           <- apply(theta_mat[, smaller_values, drop = FALSE] * binom_mat_lower, 1, max)
      lower           <- apply(theta_mat[, smaller_values, drop = FALSE], 1, max)
      
    } 
    
    # for each independent binomial proportion, the upper bound is 1
    upper <- 1
    
    xi_mat[,index[k]]  <- qnorm((theta_mat[,k] - lower)/(upper - lower))
    
  }
  
  return(xi_mat)
  
}

#' Backtransforms Samples From Real Line To Beta Parameters
#'
#' Transforms samples from the real line to samples from a truncated beta density using a stick-breaking algorithm.
#' This algorithm is suitable for mixtures of equality constrained parameters, inequality constrained parameters,
#' and free parameters
#'
#' @param xi_mat matrix with samples from truncated beta density. These samples should be transformed, so they range
#'               over the entire real line
#' @param boundaries list containing indices for upper and lower truncation boundaries
#' @param binom_equal multiplicative elements for each lower and upper bound of each inequality constrained parameter.
#' @param hyp_direction specifies whether the imposed inequality constrained imposes an increasing (i.e., 'smaller') or
#'                      decreasing (i.e., 'larger') trend
#' @return list consisting of the following elements:
#'         (1) theta_mat: matrix with transformed samples
#'         (2) lower_mat: matrix containing the lower bound for each parameter
#'         (3) upper_mat: matrix containing the upper bound for each parameter
#' @export
tbinom_backtrans      <- function(xi_mat        , boundaries, binom_equal, hyp_direction){
  # xi_mat    : samples of transformed ordered proportions (nsamples x nparam)
  # boundaries: list with upper and lower bounds
  nparameters <- length(boundaries)
  nsamples    <- nrow(xi_mat)
  theta_mat   <- matrix(NA, ncol = nparameters, nrow = nsamples)
  lower_mat   <- matrix(0 , ncol = nparameters, nrow = nsamples)
  upper_mat   <- matrix(1 , ncol = nparameters, nrow = nsamples)
  
  
  # To transform the variables we need to move from the smallest to the largest parameter.
  # Therefore, we need to determine the direction of the hypothesis and adjust the algorithm
  # accordingly.
  if(hyp_direction == 'smaller'){
    
    direction          <- 1:nparameters # move from smallest to largest parameter
    index              <- 1:nparameters # index in xi_mat matrix
    smallest_parameter <- 1
    
  } else {
    
    direction          <- nparameters:1   # move from smallest to largest parameter
    index              <- 1:nparameters   # index in xi_mat matrix
    smallest_parameter <- nparameters
    
  }
  
  for(k in direction){
    
    smaller_values <- boundaries[[k]]$lower
    larger_values  <- boundaries[[k]]$upper
    
    if(length(smaller_values) != 0){
      
      #binom_mat_lower       <- matrix(binom_equal[[k]]$lower, ncol=length(binom_equal[[k]]$lower), nrow=nrow(theta_mat), byrow=TRUE)
      #lower_mat[,index[k]] <-  apply(theta_mat[, smaller_values, drop = FALSE] * binom_mat_lower, 1, max)
      lower_mat[,index[k]] <-  apply(theta_mat[, smaller_values, drop = FALSE], 1, max)
      
    } 
    
    theta_mat[,k]        <- (upper_mat[,index[k]] - lower_mat[,index[k]]) *  pnorm(xi_mat[,index[k]]) + lower_mat[,index[k]]
    
  }
  
  return(list(theta_mat = theta_mat,
              lower_mat = lower_mat,
              upper_mat = upper_mat))
  
}