#' @title Computes Bayes Factors For Equality Constrained Multinomial Parameters
#'
#' @description Computes Bayes factor for equality constrained multinomial parameters 
#' using the standard Bayesian multinomial test.
#' Null hypothesis \eqn{H_0} states that category proportions are exactly equal to those
#' specified in \code{p}.
#' Alternative hypothesis \eqn{H_e} states that category proportions are free to vary.
#'
#' @inheritParams mult_bf_informed
#' @inherit mult_bf_informed
#' @param p numeric. A vector of probabilities of the same length as \code{x}. 
#' Its elements must be greater than 0 and less than 1. Default is 1/K
#' @return Returns a \code{data.frame} containing the Bayes factors \code{LogBFe0}, \code{BFe0}, and \code{BF0e}
#' 
#' @family functions to evaluate informed hypotheses
#' @examples 
#' data(lifestresses)
#' x <- lifestresses$stress.freq
#' a <- rep(1, nrow(lifestresses))
#' mult_bf_equality(x=x, a=a)
#' @export
mult_bf_equality <- function(x, a, p = rep(1/length(a), length(a))){
  
  # Check user input
  .checkAlphaAndData(alpha=a, counts=x)
  
  if(any(p < 0)){
    
    stop("Probabilities must be non-negative.")
    
  }
  
  if(length(x) != length(p)){
    
    stop("p and counts are not of the same length. ")
    
  }
  
  if(sum(p) != 1){
      
      p <- p/sum(p)
      warning("Parameters have been rescaled.")
      
    }
  
  # compute Bayes factor
  lbeta.xa <- sum(lgamma(a + x)) - lgamma(sum(a + x))
  lbeta.a  <- sum(lgamma(a)) - lgamma(sum(a))
  
  if (any(rowSums(cbind(p, x)) == 0)) {
    
    # in this case, x*log(p) should be zero, omit to avoid numerical issue with log(0)
    
    logBFe0 <- (lbeta.xa-lbeta.a)
    
  } else {
    
    logBFe0 <- (lbeta.xa-lbeta.a) + (0 - sum(x * log(p)))
    
  }
  
  bf <- data.frame(LogBFe0 = logBFe0,
                   BFe0    = exp(logBFe0),
                   BF0e    = 1/exp(logBFe0))
  
  return(list(bf       = bf))
  
}
