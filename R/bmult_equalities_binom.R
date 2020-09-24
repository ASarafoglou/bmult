#' Computes Bayes Factors For Equality Constrained Binomial Parameters
#'
#' Computes Bayes factor for equality constrained binomial parameters.
#' Null hypothesis H0 states that binomial proportions are exactly equal.
#' Alternative hypothesis He states that binomial proportions are free to vary.
#'
#' @inheritParams multBayesInformed
#' @inheritParams binomBayesInformed
#' @return list consisting of the following elements:
#'         (1) Bf: dataframe containing the Bayes factor logBFe0, BFe0, BF0e
#' @export
binomBfEquality <- function(a, b, counts, total){
  
  # Check user input
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