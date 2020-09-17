#' Computes Bayes Factors For Inequality Constrained Parameters
#'
#' Computes Bayes factor for inequality constrained parameters using a bridge sampling routine.
#' Restricted hypothesis Hr states that category proportions are follow a particular trend.
#' Alternative hypothesis He states that category proportions are free to vary
#'
#' @param samples martrix of dimension (nsamples x nparams) with samples from truncated Dirichlet density
#' @param restrictions either \code{character vector} containing the user specified order restriction or \code{list} of class \code{bmult_rl} as returned from \code{generateRestrictionList} that encodes 
#' inequality constraints for each independent restriction. 
#' @param index index of current restriction. Default is 1.
#' @param prior logical. If TRUE the function will ignore the data and sample from the prior distribution
#' @param maxiter maximum number of iterations for the iterative updating scheme. Default is 1,000 to avoid infinite loops.
#' @param seed set the seed for version control.
#' @param ... additional arguments (currently ignored).
#' @return list consisting of the following elements:
#'         (1) eval: list consisting of the following elements
#'                   q11: log posterior evaluations for posterior samples.
#'                   q12: log proposal evaluations for posterior samples.
#'                   q21: log posterior evaluations for samples from proposal.
#'                   q22: log proposal evaluations for samples from proposal.
#'         (2) niter: number of iterations of the iterative updating scheme.
#'         (3) logml: estimate of log marginal likelihood.
#'         (4) hyp: character vector that contains the inequality constrained hypothesis 
#' @export
multBayesBfInequality <- function(samples, restrictions, alpha = rep(1,ncol(samples)), data = NULL, prior = FALSE, 
                                  index = 1, maxiter = 1e3, seed=NULL, ...){
  
  ###    Code by Gronau et al. (2017) - online appendix ###
  ###    Modified by Alexandra Sarafoglou               ###
  
  # Note that before applying this function the user needs to:
  # 1. Collect 2*N1 samples from the truncated prior and posterior distribution
  #    (e.g., through MCMC sampling)
  # 2. Choose a suitable proposal distribution. Here we choose the multivariate normal &
  #    Specify the function for evaluating the log of the unnormalized density.
  #    This function is here referred to as log_unnormalized_density.
  
  # 0.1 Check User Input
  .checksIfMatrix(samples)
  
  ## 0.2 Extract Relevant Information
  .checkRestrictionListClass(restrictions)
  
  # if order restriction is given as character vector, create restriction list
  if(!inherits(restrictions, 'bmult_rl') & !inherits(restrictions, 'bmult_rl_ineq')){
    
    if(!is.null(colnames(samples))){
      
      factor_levels <- colnames(samples)
      
    } else {
      
      factor_levels <- paste0('theta', 1:ncol(samples))
      
    }
    
    restriction_list <- generateRestrictionList(restrictions, factor_levels, alpha, data)
    # only consider inequality constraints
    restrictions     <- restriction_list$inequality_constraints
    prior_and_data   <- restrictions$alpha_inequalities[[index]] + restrictions$counts_inequalities[[index]]
    
  } else {
    
    if(inherits(restrictions, 'bmult_rl')){
      
      restrictions <- restrictions$inequality_constraints
      
    }
    
    if(prior){
      
      prior_and_data   <- restrictions$alpha_inequalities[[index]]
      
    } else {
      
      prior_and_data   <- restrictions$alpha_inequalities[[index]] + restrictions$counts_inequalities[[index]]
      
    }
    
  }
  
  boundaries       <- restrictions$boundaries[[index]]
  nr_mult_free     <- restrictions$nr_mult_free[[index]]
  nr_mult_equal    <- restrictions$nr_mult_equal[[index]]
  mult_equal       <- restrictions$mult_equal[[index]]
  hyp_direction    <- restrictions$direction[index]
  hyp              <- restrictions$hyp[[index]]
  
  # set seed if wanted
  if(!is.null(seed) & is.numeric(seed)){
    set.seed(seed)
  }
  
  # check if correct number of parameters were provided
  .checkNrParameters(samples = samples, boundaries = boundaries, data = prior_and_data)
  
  # 2. Specify the function for evaluating the log of the unnormalized density
  # 3. Transform the parameters to the real line
  samples     <- tDirTrans(samples, boundaries, mult_equal, nr_mult_equal, nr_mult_free, hyp_direction)
  # 4. Split the samples into two parts
  # Use the first 50% for fiting the proposal distribution and the second 50%
  # in the iterative scheme.
  nperchain      <- nrow(samples)
  fit_index      <- 1:(nperchain/2)
  samples_4_fit  <- samples[fit_index,, drop = FALSE]
  samples_4_iter <- samples[-fit_index,, drop = FALSE]
  
  # 5. Fit proposal distribution
  N2 <- N1 <- nrow(samples_4_iter)
  m  <- apply(samples_4_fit, 2, mean) # mean vector
  V  <- cov(samples_4_fit)            # covariance matrix
  # 6. Draw N2 samples from the proposal distribution
  gen_samples <- mvtnorm::rmvnorm(N2, m, V)
  # 7a. Evaluate proposal distribution for posterior & generated samples
  q12 <- mvtnorm::dmvnorm(samples_4_iter, m, V, log = TRUE)
  q22 <- mvtnorm::dmvnorm(gen_samples   , m, V, log = TRUE)
  # 7b. Evaluate unnormalized posterior for posterior & generated samples
  q11 <- logUnnormalizedTDir(samples_4_iter, boundaries, mult_equal, nr_mult_equal, nr_mult_free, hyp_direction, prior_and_data)
  q21 <- logUnnormalizedTDir(gen_samples   , boundaries, mult_equal, nr_mult_equal, nr_mult_free, hyp_direction, prior_and_data)
  
  # 8. Run iterative scheme as proposed in Meng and Wong (1996) to estimate
  # the marginal likelihood
  l1 <- q11 - q12
  l2 <- q21 - q22
  # increase numerical stability by subtracting the median of l1 from l1 & l2
  lstar <- median(l1)
  s1    <- N1/(N1 + N2)
  s2    <- N2/(N1 + N2)
  e     <- Brobdingnag::as.brob( exp(1) )     # more stable Brobdingnag number representation
  criterion_val <- 1e-10 + 1 # criterion value
  r <- 0                     # starting value for r
  i <- 0                     # iteration counter
  
  while (criterion_val > 1e-10 & i < maxiter) {
    r_old <- r
    numerator <- as.numeric(e^(l2 - lstar)/(s1 * e^(l2 - lstar) + s2 *  r))
    denominator <- as.numeric(1/(s1 * e^(l1 - lstar) + s2 * r))
    r <- (N1/N2)*sum(numerator)/sum(denominator)
    i <- i + 1
    criterion_val <- abs((r - r_old)/r)
  }
  
  logml <- log(r) + lstar # log of marginal likelihood
  
  # Return a list with the evaluations of the proposal and the unnormalized
  # posterior, the number of iterations of the iterative scheme, and the
  # estimated log marginal likelihood
  output <- list(eval  = list(q11 = q11, q12 = q12,
                              q21 = q21, q22 = q22),
                 niter = i, logml = logml, hyp = hyp)
  # Compute error measures for estimated marginal likelihood
  error_measures        <- .computeRMSE(output)
  output$error_measures <- error_measures
  
  # assign class
  class(output) <- 'bmult_bridge'
  return(output)
}