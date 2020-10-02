#' Computes The Relative Mean-Squared Error For Marginal Likelihood Estimate
#'
#' Computes error measures for estimated marginal likelihood.
#'
#' @usage 
#' \code{
#' computeRMSE(bridge_object, ...)
#' ## S3 method for class 'bridge'
#' computeRMSE(bridge_object, ...)
#' }
#' @param bridge_output an object of class \code{"bridge"} as returned from \code{multBfInformed}
#' @param ... additional arguments (currently ignored).
#' @details Computes error measures for marginal likelihood bridge sampling estimates. The approximate errors for bridge sampling estimates are based on Fruehwirth-Schnatter (2004). 
#' Code is based on \code{error_measures} function of the \code{R} package \code{bridgesampling}.
#' @author Quentin F. Gronau
#' @references \insertRef{gronau2017bridgesampling}{multibridge} \insertRef{fruhwirth2004estimating}{multibridge}
#' @return Returns a list with the following components:
#' \describe{
#'   \item{re2}{approximate relative mean-squared error for marginal likelihood estimate.}
#'   \item{cv}{approximate coefficient of variation for marginal likelihood estimate (assumes that bridge estimate is unbiased).}
#'   \item{percentage}{approximate percentage error of marginal likelihood estimate.}
#' }
.computeRMSE <- function(bridge_output) {
  
  # function that computes an approximate relative mean-squared error for
  # a marginal likelihood estimated via bridge sampling
  # (see Fruehwirth-Schnatter, 2004)
  # Code by Gronau et al. (2017)
  
  e <- Brobdingnag::as.brob( exp(1) )
  
  ml <- e^(bridge_output$logml)
  g_p <- e^(bridge_output$eval$q12)
  g_g <- e^(bridge_output$eval$q22)
  priorTimesLik_p <- e^(bridge_output$eval$q11)
  priorTimesLik_g <- e^(bridge_output$eval$q21)
  p_p <- priorTimesLik_p/ml
  p_g <- priorTimesLik_g/ml
  
  N1 <- length(p_p)
  N2 <- length(g_g)
  s1 <- N1/(N1 + N2)
  s2 <- N2/(N1 + N2)
  
  f1 <- as.numeric( p_g/(s1*p_g + s2*g_g) )
  f2 <- as.numeric( g_p/(s1*p_p + s2*g_p) )
  rho_f2 <- coda::spectrum0.ar( f2 )$spec
  
  term1 <- 1/N2 * var( f1 ) / mean( f1 )^2
  term2 <- rho_f2/N1 * var( f2 ) / mean( f2 )^2
  
  re2        <- term1 + term2
  cv         <- sqrt(re2)
  percentage <- paste0(round(cv*100, 4), '%')
  
  output <- list(re2=re2,
                 cv=cv,
                 percentage=percentage)
  return(output)
  
}
