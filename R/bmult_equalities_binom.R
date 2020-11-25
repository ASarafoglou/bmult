#' @title Computes Bayes Factors For Equality Constrained Binomial Parameters
#'
#' @description Computes Bayes factor for equality constrained binomial parameters.
#' Null hypothesis \eqn{H_0} states that binomial proportions are exactly equal.
#' Alternative hypothesis \eqn{H_e} states that binomial proportions are free to vary.
#'
#' @inherit binom_bf_informed
#' @inheritParams binom_bf_informed
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
binom_bf_equality <- function(x, n=NULL, a, b){
  
  # Check user input
  userInput <- .checkIfXIsVectorOrTable(x, n)
  counts    <- userInput$counts
  total     <- userInput$total
  .checkAlphaAndData(alpha=a, beta=b, counts=counts, total=total)
  
  # compute Bayes factor
  
  lbeta.xa.H0 <- lbeta(sum(counts) + sum(a) - length(a) + 1, sum(total) - sum(counts) + sum(b) - length(b) + 1)
  lbeta.a.H0  <- lbeta(sum(a) - length(a) + 1, sum(b) - length(b) + 1)
  
  lbeta.xa.He <- sum(lbeta(counts + a, total - counts + b ))
  lbeta.a.He  <- sum(lbeta(a, b))
  
  logBFe0 <-  (lbeta.xa.He-lbeta.a.He) - (lbeta.xa.H0-lbeta.a.H0)
  
  bf <- data.frame(LogBFe0 = logBFe0,
                   BFe0    = exp(logBFe0),
                   BF0e    = 1/exp(logBFe0))
  
  return(list(bf = bf))
  
}