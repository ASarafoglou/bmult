#' Computes Length Of Remaining Stick
#'
#' When applying the probit transformation on the Dirichlet samples, this function is used as part of the stick-breaking
#' algorithm. It computes the length of the remaining stick when the current element is broken off.
#'
#' @param theta_mat matrix with samples from truncated Dirichlet density
#' @param k current parameter index
#' @param hyp_direction specifies whether the imposed inequality constrained imposes an increasing (i.e., 'smaller') or
#'                      decreasing (i.e., 'larger') trend
#' @return stick length
.computeLengthOfRemainingStick       <- function(theta_mat     , k         , hyp_direction) {
  if(hyp_direction == 'smaller'){
    stick_length   <- (1 - rowSums(theta_mat[, 1:(k - 1), drop = FALSE]))
  } else {
    K              <- ncol(theta_mat)
    stick_length   <- (1 - rowSums(theta_mat[, (k + 1):K, drop = FALSE]))
  }
  return(stick_length)
}

#' Adjusts Upper Bound For Free Parameters
#'
#' Corrects the upper bound for current parameter. This correction only applies for parameters that are free to vary within
#' the restriction. Then the length of the remaining stick must be based on the largest free parameter value.
#'
#' @param theta_mat matrix with samples from truncated Dirichlet density
#' @param k current parameter index
#' @param upper current upper bound
#' @param nr_mult_equal vector of multiplicative elements of collapsed parameters
#' @param smaller_values index of parameters that are smaller than the current one
#' @param larger_values index of parameters that are larger than the current one
#' @param hyp_direction specifies whether the imposed inequality constrained imposes an increasing (i.e., 'smaller') or
#'                      decreasing (i.e., 'larger') trend
#' @return adjusted upper bound
.adjustUpperBoundForFreeParameters   <- function(theta_mat     , k, upper, nr_mult_equal, smaller_values, larger_values, hyp_direction){
  
  # only select free parameters that have been evaluated already
  free_parameters   <- seq_along(nr_mult_equal)[-c(smaller_values, larger_values)]
  
  if(hyp_direction == 'smaller'){
    
    smaller_parameters <- 1:(k - 1)
    
  } else {
    
    K                  <- length(nr_mult_equal)
    smaller_parameters <- k:K
    
  }
  
  free_parameters   <- free_parameters[free_parameters %in% smaller_parameters]
  delete_NA_columns <- apply(theta_mat[, free_parameters, drop = FALSE], 2, function(x) any(is.na(x)))
  free_parameters   <- free_parameters[!delete_NA_columns]
  
  if(length(free_parameters) != 0) {
    
    # remove the correction for equality constraints to get the real boundaries
    mult_mat_free_raw  <- matrix(1/nr_mult_equal[free_parameters], ncol=length(nr_mult_equal[free_parameters]), nrow=nrow(theta_mat), byrow=TRUE)
    max_free_raw       <- apply(theta_mat[, free_parameters, drop = FALSE] * mult_mat_free_raw, 1, max)
    adjustment         <- ifelse(max_free_raw > upper, max_free_raw - upper, 0)
    adjustment         <- adjustment * sum(nr_mult_equal[larger_values]) * (1/nr_mult_equal[k])
    upper              <- upper - adjustment
    
  }
  
  return(upper)
  
}

#' Transforms Truncated Dirichlet Samples To Real Line
#'
#' Transforms samples from a truncated Dirichlet density to the real line using a stick-breaking algorithm.
#' This algorithm is suitable for mixtures of equality constrained parameters, inequality constrained parameters,
#' and free parameters
#'
#' @param theta_mat matrix with samples from truncated Dirichlet density
#' @param boundaries list containing indices for upper and lower truncation boundaries
#' @param mult_equal multiplicative elements for each lower and upper bound of each inequality constrained parameter.
#' @param nr_mult_equal vector of multiplicative elements of collapsed parameters
#' @param nr_mult_free vector of multiplicative elements of free parameters
#' @param hyp_direction specifies whether the imposed inequality constrained imposes an increasing (i.e., 'smaller') or
#'                      decreasing (i.e., 'larger') trend
#' @return matrix with transformed samples
#' @export
tdir_trans           <- function(theta_mat     , boundaries, mult_equal, nr_mult_equal, nr_mult_free, hyp_direction){
  # theta_mat : prior or posterior samples (nsamples x nparam)
  # boundaries: list with upper and lower bounds
  nparameters   <- length(boundaries)
  nsamples      <- nrow(theta_mat)
  xi_mat        <- matrix(NA, ncol = c(nparameters - 1), nrow = nsamples)
  
  # To transform the variables we need to move from the smallest to the largest parameter.
  # Therefore, we need to determine the direction of the hypothesis and adjust the algorithm
  # accordingly.
  if(hyp_direction == 'smaller'){
    
    direction          <- 1:(nparameters - 1) # move from smallest to largest parameter
    index              <- 1:(nparameters - 1) # index in xi_mat matrix
    smallest_parameter <- 1
    
  } else {
    
    direction          <- (nparameters):2           # move from smallest to largest parameter
    index              <- c(NA,1:(nparameters-1))   # index in xi_mat matrix
    smallest_parameter <- nparameters
    
  }
  
  for(k in direction){
    
    smaller_values <- boundaries[[k]]$lower
    larger_values  <- boundaries[[k]]$upper
    
    if(length(smaller_values) != 0){
      
      mult_mat_lower   <- matrix(mult_equal[[k]]$lower, ncol=length(mult_equal[[k]]$lower), nrow=nrow(theta_mat), byrow=TRUE)
      lower            <- apply(theta_mat[, smaller_values, drop = FALSE] * mult_mat_lower, 1, max)
      
      # adjust for following parameter
      mult_mat_lower_raw <- matrix(1/nr_mult_equal[smaller_values], ncol=length(nr_mult_equal[smaller_values]), nrow=nrow(theta_mat), byrow=TRUE)
      lower_raw          <- apply(theta_mat[, smaller_values, drop = FALSE] * mult_mat_lower_raw, 1, max)
      
    } else {
      
      lower_raw <- rep(0, nsamples)
      lower     <- rep(0, nsamples)
      
    }
    
    if(k == smallest_parameter){
      
      length_remaining_stick   <- 1
      
    } else {
      
      length_remaining_stick <- .computeLengthOfRemainingStick(theta_mat, k, hyp_direction)
      # adjust for free parameters; subtract lower bound of remaining free parameters from stick length
      length_remaining_stick <- length_remaining_stick - (nr_mult_free[k] * lower_raw)
      
    }
    
    # elements remaining stick: nr_larger values + current value
    elements_remaining_stick <- 1 * nr_mult_equal[k] + sum(nr_mult_equal[larger_values])
    
    upper <- (length_remaining_stick/elements_remaining_stick)
    # adjust for free parameters; if free parameters are in the constraint, and they are larger
    # than the current one, adjust for the discrepancy
    upper <- .adjustUpperBoundForFreeParameters(theta_mat, k, upper, nr_mult_equal, smaller_values, larger_values, hyp_direction)
    
    # adjust for equal parameters; if current parameter is equality constrained, enlarge upper bound
    upper <- upper * nr_mult_equal[k]
    
    xi_mat[,index[k]]  <- qnorm((theta_mat[,k] - lower)/(upper - lower))
    
  }
  
  return(xi_mat)
  
}

#' Backtransforms Samples From Real Line To Dirichlet Parameters
#'
#' Transforms samples from the real line to samples from a truncated Dirichlet density using a stick-breaking algorithm.
#' This algorithm is suitable for mixtures of equality constrained parameters, inequality constrained parameters,
#' and free parameters
#'
#' @param xi_mat matrix with samples from truncated Dirichlet density. These samples should be transformed, so they range
#'               over the entire real line
#' @param boundaries list containing indices for upper and lower truncation boundaries
#' @param mult_equal multiplicative elements for each lower and upper bound of each inequality constrained parameter.
#' @param nr_mult_equal vector of multiplicative elements of collapsed parameters
#' @param nr_mult_free vector of multiplicative elements of free parameters
#' @param hyp_direction specifies whether the imposed inequality constrained imposes an increasing (i.e., 'smaller') or
#'                      decreasing (i.e., 'larger') trend
#' @return list consisting of the following elements:
#'         (1) theta_mat: matrix with transformed samples
#'         (2) lower_mat: matrix containing the lower bound for each parameter
#'         (3) upper_mat: matrix containing the upper bound for each parameter
#' @export
tdir_backtrans      <- function(xi_mat        , boundaries, mult_equal, nr_mult_equal, nr_mult_free, hyp_direction){
  # xi_mat    : samples of transformed ordered proportions (nsamples x nparam)
  # boundaries: list with upper and lower bounds
  nparameters            <- length(boundaries)
  nsamples               <- nrow(xi_mat)
  theta_mat              <- matrix(NA, ncol = nparameters    , nrow = nsamples)
  lower_mat <- upper_mat <- matrix(NA, ncol = nparameters - 1, nrow = nsamples)
  lower_mat_raw          <- matrix(NA, ncol = nparameters - 1, nrow = nsamples)
  
  # To transform the variables we need to move from the smallest to the largest parameter.
  # Therefore, we need to determine the direction of the hypothesis and adjust the algorithm
  # accordingly.
  if(hyp_direction == 'smaller'){
    
    direction          <- 1:(nparameters - 1) # move from smallest to largest parameter
    index              <- 1:(nparameters - 1) # index in xi_mat matrix
    smallest_parameter <- 1
    largest_parameter  <- nparameters
    
  } else {
    
    direction          <- (nparameters):2           # move from smallest to largest parameter
    index              <- c(NA,1:(nparameters-1))   # index in xi_mat matrix
    smallest_parameter <- nparameters
    largest_parameter  <- 1
  }
  
  for(k in direction){
    
    smaller_values <- boundaries[[k]]$lower
    larger_values  <- boundaries[[k]]$upper
    
    if(length(smaller_values) != 0){
      
      mult_mat_lower       <- matrix(mult_equal[[k]]$lower, ncol=length(mult_equal[[k]]$lower), nrow=nrow(theta_mat), byrow=TRUE)
      lower_mat[,index[k]] <-  apply(theta_mat[, smaller_values, drop = FALSE] * mult_mat_lower, 1, max)
      
      # adjust for following parameter
      mult_mat_lower_raw       <- matrix(1/nr_mult_equal[smaller_values], ncol=length(nr_mult_equal[smaller_values]), nrow=nrow(theta_mat), byrow=TRUE)
      lower_mat_raw[,index[k]] <- apply(theta_mat[, smaller_values, drop = FALSE] * mult_mat_lower_raw, 1, max)
      
    } else {
      
      lower_mat_raw[,index[k]] <- rep(0, nsamples)
      lower_mat[,index[k]]     <- rep(0, nsamples)
      
    }
    
    # elements remaining stick: nr_larger values + current value
    elements_remaining_stick <- 1 * nr_mult_equal[k] + sum(nr_mult_equal[larger_values])
    
    if(k == smallest_parameter){
      
      length_remaining_stick <- 1
      
    } else {
      
      length_remaining_stick <- .computeLengthOfRemainingStick(theta_mat, k, hyp_direction)
      # adjust for free parameters; subtract previous treated free parameters from stick length
      length_remaining_stick <- length_remaining_stick - (nr_mult_free[k] * lower_mat_raw[,index[k]])
      
    }
    
    upper_mat[,index[k]] <- (length_remaining_stick/elements_remaining_stick)
    # adjust for free parameters; if free parameters are in the constraint, and they are larger
    # than the current one, adjust for the discrepancy
    upper_mat[,index[k]] <- .adjustUpperBoundForFreeParameters(theta_mat, k, upper_mat[,index[k]], nr_mult_equal, smaller_values, larger_values, hyp_direction)
    
    # adjust for equal parameters; if current parameter is equality constrained, enlarge upper bound
    upper_mat[,index[k]] <- upper_mat[,index[k]] * nr_mult_equal[k]
    
    theta_mat[,k]        <- (upper_mat[,index[k]] - lower_mat[,index[k]]) *  pnorm(xi_mat[,index[k]]) + lower_mat[,index[k]]
    
  }
  # sum-to-one constraint applies to lambda, not lambda/2
  # Therefore, multiply theta with the number of multiplicative elements
  theta_mat[, largest_parameter] <- 1 - rowSums(theta_mat[,, drop = FALSE], na.rm = TRUE)
  
  return(list(theta_mat = theta_mat,
              lower_mat = lower_mat,
              upper_mat = upper_mat))
  
}