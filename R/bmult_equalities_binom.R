#' @title Computes Bayes Factors For Equality Constrained Binomial Parameters
#'
#' @description Computes Bayes factor for equality constrained binomial parameters.
#' Null hypothesis \eqn{H_0} states that binomial proportions are exactly equal or
#' exactly equal and equal to \code{p}.
#' Alternative hypothesis \eqn{H_e} states that binomial proportions are free to vary.
#'
#' @inherit binom_bf_informed
#' @inheritParams binom_bf_informed
#' @param p numeric. Hypothesized probability of success. Must be greater than 0 and less than 1.
#' Default sets all binomial proportions exactly equal without specifying a specific value.
#' @return Returns a \code{data.frame} containing the Bayes factors \code{LogBFe0}, \code{BFe0}, and \code{BF0e}
#' 
#' @family functions to evaluate informed hypotheses
#' @examples 
#' data(journals)
#' x <- journals$errors
#' n <- journals$nr_NHST
#' a <- rep(1, nrow(journals))
#' b <- rep(1, nrow(journals))
#' binom_bf_equality(x=x, n=n, a=a, b=b)
#' @export
binom_bf_equality <- function(x, n=NULL, a, b, p = NULL){
  
  # Check user input
  counts    <- x
  total     <- n
  
  # compute Bayes factor
  lbeta.xa.He <- sum(lbeta(counts + a, total - counts + b ))
  lbeta.a.He  <- sum(lbeta(a, b))
  
  # all underlying probabilities are equal
  if(is.null(p)){
    
    lbeta.xa.H0 <- lbeta(sum(counts) + sum(a) - length(a) + 1, sum(total) - sum(counts) + sum(b) - length(b) + 1)
    lbeta.a.H0  <- lbeta(sum(a) - length(a) + 1, sum(b) - length(b) + 1)
    logBFe0     <-  (lbeta.xa.He-lbeta.a.He) - (lbeta.xa.H0-lbeta.a.H0)
    
  } else {
    # all underlying probabilities are equal and equal to specific values
    
    if (p == 0 && sum(counts) == 0) {
      
      # in this case, counts*log(p) should be zero, omit to avoid numerical issue with log(0)
      loglikelihood.H0 <- log(1 - p)*(sum(total) - sum(counts))
      
    } else if (p == 1 && sum(counts) == sum(total)) {
      
      # in this case, (n - counts)*log(1 - p) should be zero, omit to avoid numerical issue with log(0)
      loglikelihood.H0 <- log(p)*sum(counts)
      
    } else {
      
      loglikelihood.H0 <- log(p)*sum(counts) + log(1 - p)*(sum(total) - sum(counts))
    }
    
    logBFe0 <- (lbeta.xa.He - lbeta.a.He) - loglikelihood.H0
    
  }
  
  bf <- data.frame(LogBFe0 = logBFe0,
                   BFe0    = exp(logBFe0),
                   BF0e    = 1/exp(logBFe0))
  
  return(list(bf = bf))
  
}
