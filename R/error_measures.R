# Computes The Relative Mean-Squared Error For Marginal Likelihood Estimate
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
  
  term1 <- 1/N2 * stats::var( f1 ) / mean( f1 )^2
  term2 <- rho_f2/N1 * stats::var( f2 ) / mean( f2 )^2
  
  re2        <- term1 + term2
  cv         <- sqrt(re2)
  percentage <- paste0(round(cv*100, 4), '%')
  
  output <- list(re2=re2,
                 cv=cv,
                 percentage=percentage)
  return(output)
  
}
