#' Computes Bayes Factors For Inequality Constrained Inedepent Binomial Parameters
#'
#' Computes Bayes factor for inequality constrained binomial parameters using a bridge sampling routine.
#' Restricted hypothesis Hr states that binomial proportions follow a particular trend.
#' Alternative hypothesis He states that binomial proportions are free to vary.
#'
#' @param samples matrix of dimension (nsamples x nparams) with samples from independent truncated Binomial densities
#' @param restrictions either \code{character vector} containing the user specified order restriction or \code{list} of class \code{bmult_rl} as returned from \code{generateRestrictionList} that encodes 
#' inequality constraints for each independent restriction.
#' @param a numeric vector with values of alpha parameters of the beta distribution 
#' @param b numeric vector with values of beta parameters of the beta distribution 
#' @param counts numeric vector with number of successes
#' @param total numeric vector with total number of observations
#' @param prior logical. If TRUE the function will ignore the data and sample from the prior distribution
#' @param index index of current restriction. Default is 1.
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
binomBfInequality <- function(samples, restrictions, a = rep(1,ncol(samples)), b = rep(1,ncol(samples)), counts = NULL, total = NULL, prior = FALSE, 
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
    
    .checkAlphaAndData(alpha=a, beta=b, counts=counts, total=total)
    restriction_list <- generateRestrictionList(restrictions, factor_levels, a=a, b=b, counts=counts, total=total, binom=TRUE)
    # only consider inequality constraints
    restrictions     <- restriction_list$inequality_constraints
    a                <- restrictions$alpha_inequalities[[index]] + restrictions$counts_inequalities[[index]]
    b                <- restrictions$beta_inequalities[[index]] + (restrictions$total_inequalities[[index]] - restrictions$counts_inequalities[[index]])
    
  } else {
    
    if(inherits(restrictions, 'bmult_rl')){
      
      restrictions <- restrictions$inequality_constraints
      
    }
    
    if(prior){
      
      a   <- restrictions$alpha_inequalities[[index]]
      b   <- restrictions$beta_inequalities[[index]]
      
    } else {
      
      a   <- restrictions$alpha_inequalities[[index]] + restrictions$counts_inequalities[[index]]
      b   <- restrictions$beta_inequalities[[index]] + (restrictions$total_inequalities[[index]] - restrictions$counts_inequalities[[index]])
      
    }
    
  }
  
  boundaries       <- restrictions$boundaries[[index]]
  binom_equal      <- restrictions$mult_equal[[index]]
  hyp_direction    <- restrictions$direction[index]
  hyp              <- restrictions$hyp[[index]]
  
  # set seed if wanted
  if(!is.null(seed) & is.numeric(seed)){
    set.seed(seed)
  }
  
  # check if correct number of parameters were provided
  .checkNrParameters(samples = samples, boundaries = boundaries)
  
  # 2. Specify the function for evaluating the log of the unnormalized density
  # 3. Transform the parameters to the real line
  samples     <- tBinomTrans(samples, boundaries, binom_equal, hyp_direction)
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
  q11 <- logUnnormalizedTBinom(samples_4_iter, boundaries, binom_equal, hyp_direction, a, b)
  q21 <- logUnnormalizedTBinom(gen_samples   , boundaries, binom_equal, hyp_direction, a, b)
  
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

#' #' Computes Bayes Factors For Inequality Constrained Binomial Parameters
#' #'
#' #' Computes Bayes factor for inequality constrained binomial parameters using the encompassing prior approach.
#' #' Restricted hypothesis Hr states that binomial proportions are follow a particular trend.
#' #' Alternative hypothesis He states that binomial proportions are free to vary.
#' #'
#' #' @param samples matrix of dimension (nsamples x nparams) with samples from encompassing binomial densities
#' #' @param restrictions either \code{character vector} containing the user specified order restriction or \code{list} of class \code{bmult_rl} as returned from \code{generateRestrictionList} that encodes 
#' #' inequality constraints for each independent restriction. 
#' #' @param a numeric vector with values of alpha parameters of the beta distribution 
#' #' @param b numeric vector with values of beta parameters of the beta distribution 
#' #' @param counts numeric vector with number of successes
#' #' @param total numeric vector with total number of observations
#' #' @param index index of current restriction. Default is 1.
#' #' @param ... additional arguments (currently ignored).
#' #' @return proportion of samples in accordance with the constraint (in log scale for numeric stability)
#' #' @export
#' binomBfInequality    <- function(samples, restrictions, a = rep(1,ncol(samples)), b = rep(1,ncol(samples)), counts = NULL, total = NULL, index  = 1, ...){
#'   
#'   # 0.1 Check User Input
#'   .checksIfMatrix(samples)
#'   
#'   ## 0.2 Extract Relevant Information
#'   .checkRestrictionListClass(restrictions)
#'   
#'   # if order restriction is given as character vector, create restriction lis 
#'   # if order restriction is given as character vector, create restriction list
#'   if(!inherits(restrictions, 'bmult_rl') & !inherits(restrictions, 'bmult_rl_ineq')){
#'     
#'     if(!is.null(colnames(samples))){
#'       
#'       factor_levels <- colnames(samples)
#'       
#'     } else {
#'       
#'       factor_levels <- paste0('theta', 1:ncol(samples))
#'       
#'     }
#'     
#'     # before creating restriction list
#'     .checkAlphaAndData(alpha=a, beta=b, counts=counts, total=total)
#'     restriction_list <- generateRestrictionList(restrictions, factor_levels, a=a, b=b, counts=counts, total=total, binom=TRUE)
#'     # only consider inequality constraints
#'     restrictions     <- restriction_list$inequality_constraints
#'     
#'   } else {
#'     
#'     if(inherits(restrictions, 'bmult_rl')){
#'       
#'       restrictions <- restrictions$inequality_constraints
#'       
#'     }
#'     
#'   }
#'   
#'   boundaries       <- restrictions$boundaries[[index]]
#'   
#'   # check if correct number of parameters were provided
#'   .checkNrParameters(samples = samples, boundaries = boundaries)
#'   
#'   N        <- ncol(samples)
#'   nsamples <- nrow(samples)
#'   I_mat    <- matrix(NA, nrow = nsamples, ncol = N)
#'   
#'   # initialize progress bar
#'   progress_bar_text  <- paste0("restr. ", index , ". evaluation completed: [:bar] :percent time remaining: :eta")
#'   pb  <- progress::progress_bar$new(format = progress_bar_text, total = 2*N, clear = FALSE, width= 80)
#'   
#'   for(i in 1:N){
#'     
#'     obeys_lower <- rowSums(samples[, i] < samples[, boundaries[[i]]$lower,  drop = FALSE]) > 0
#'     pb$tick()
#'     obeys_upper <- rowSums(samples[, i] > samples[, boundaries[[i]]$upper,  drop = FALSE]) > 0
#'     pb$tick()
#'     
#'     I_mat[, i]  <- obeys_lower & obeys_upper
#'     
#'   }
#'   I <- mean(rowSums(I_mat) == N)
#'   
#'   # Returns the (log) proportion of samples in accordance with the constraint
#'   output <- log(I)
#'   
#'   return(output)
#'   
#' }