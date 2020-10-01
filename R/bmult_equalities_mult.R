#' @title Computes Bayes Factors For Equality Constrained Multinomial Parameters
#'
#' @description Computes Bayes factor for equality constrained multinomial parameters 
#' using the standard Bayesian multinomial test.
#' Null hypothesis \eqn{H_0} states that category proportions are exactly equal to predefined values.
#' Alternative hypothesis \eqn{H_e} states that category proportions are free to vary.
#'
#' @inheritParams multBayesInformed
#' @param theta numeric. Values of interest. Default is 1/K
#' @return list consisting of the following elements:
#' \describe{
#' \item{\code{$bf}}{\code{data.frame} containing the Bayes factors \code{LogBFe0}, \code{BFe0}, and \code{BF0e}}
#' \item{\code{$expected}}{numeric. vector with expected values}
#' }
#' @family functions to evaluate informed hypotheses
#' @examples 
#' data(lifestresses)
#' x <- lifestresses$stress.freq
#' a <- rep(1, nrow(lifestresses))
#' multBfEquality(x=x, a=a)
#' @export
multBfEquality <- function(x, a, theta = rep(1/length(a), length(a))){
  
  # Check user input
  .checkAlphaAndData(alpha=a, counts=x)
    
    if(sum(theta) != 1){
      theta <- theta/sum(theta)
      warning("Parameters have been rescaled.")
    }
  
  expected <- sum(x)*theta
  # compute Bayes factor
  lbeta.xa <- sum(lgamma(a + x)) - lgamma(sum(a + x))
  lbeta.a  <- sum(lgamma(a)) - lgamma(sum(a))
  
  if (any(rowSums(cbind(theta, x)) == 0)) {
    
    # in this case, x*log(theta) should be zero, omit to avoid numerical issue with log(0)
    
    logBFe0 <- (lbeta.xa-lbeta.a)
    
  } else {
    
    logBFe0 <- (lbeta.xa-lbeta.a) + (0 - sum(x * log(theta)))
    
  }
  
  bf <- data.frame(LogBFe0 = logBFe0,
                   BFe0    = exp(logBFe0),
                   BF0e    = 1/exp(logBFe0))
  
  return(list(bf       = bf,
              expected = expected))
  
}
