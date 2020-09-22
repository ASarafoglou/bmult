#' Computes Bayes Factors For Equality Constrained Parameters
#'
#' Computes Bayes factor for equality constrained parameters using the standard Bayesian multinomial test.
#' Null hypothesis H0 states that category proportions are exactly equal to predefined values.
#' Alternative hypothesis He states that category proportions are free to vary
#'
#' @param a vector with concentration parameters
#' @param counts vector with data
#' @param theta values of interest. Default is 1/K
#' @param ... additional arguments (currently ignored).
#' @return list consisting of the following elements:
#'         (1) Bf: dataframe containing the Bayes factor logBFe0, BFe0,
#'         (3) expected: vector with expected values
#' @export
multBfEquality <- function(a, counts, theta = rep(1/length(a), length(a)), ...){
  
  # Check user input
  .checkAlphaAndData(a=a, counts=counts)
    
    if(sum(theta) != 1){
      theta <- theta/sum(theta)
      warning("Parameters have been rescaled.")
    }
  
  expected <- sum(counts)*theta
  # compute Bayes factor
  lbeta.xa <- sum(lgamma(a + counts)) - lgamma(sum(a + counts))
  lbeta.a  <- sum(lgamma(a)) - lgamma(sum(a))
  
  if (any(rowSums(cbind(theta, counts)) == 0)) {
    
    # in this case, counts*log(theta) should be zero, omit to avoid numerical issue with log(0)
    
    logBFe0 <- (lbeta.xa-lbeta.a)
    
  } else {
    
    logBFe0 <- (lbeta.xa-lbeta.a) + (0 - sum(counts * log(theta)))
    
  }
  
  bf <- data.frame(LogBFe0 = logBFe0,
                   BFe0    = exp(logBFe0),
                   BF0e    = 1/exp(logBFe0))
  
  return(list(bf       = bf,
              expected = expected))
  
}