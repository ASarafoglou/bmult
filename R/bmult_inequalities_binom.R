#' @title Computes Bayes Factors For Inequality Constrained Independent Binomial Parameters
#'
#' @description Computes Bayes factor for inequality constrained binomial parameters using a bridge sampling routine.
#' Restricted hypothesis \eqn{H_r} states that binomial proportions follow a particular trend.
#' Alternative hypothesis \eqn{H_e} states that binomial proportions are free to vary.
#' 
#' @inherit binom_bf_informed
#' @inherit mult_bf_inequality
#' @inheritParams binom_bf_informed
#' @inheritParams mult_bf_informed
#' @param samples matrix of dimension (\code{nsamples x nparams}) with samples from truncated Dirichlet density
#' @note 
#' The following signs can be used to encode restricted hypotheses: \code{"<"} and \code{">"} for inequality constraints, \code{"="} for equality constraints,
#' \code{","} for free parameters, and \code{"&"} for independent hypotheses. The restricted hypothesis can either be a string or a character vector.
#' For instance, the hypothesis \code{c("theta1 < theta2, theta3")} means 
#' \itemize{
#' \item \code{theta1} is smaller than both \code{theta2} and \code{theta3}
#' \item The parameters \code{theta2} and \code{theta3} both have \code{theta1} as lower bound, but are not influenced by each other.
#' }
#' The hypothesis \code{c("theta1 < theta2 = theta3 & theta4 > theta5")} means that 
#' \itemize{
#' \item Two independent hypotheses are stipulated: \code{"theta1 < theta2 = theta3"} and \code{"theta4 > theta5"}
#' \item The restrictions on the parameters \code{theta1}, \code{theta2}, and \code{theta3} do
#' not influence the restrictions on the parameters \code{theta4} and \code{theta5}.
#' \item \code{theta1} is smaller than \code{theta2} and \code{theta3}
#' \item \code{theta2} and \code{theta3} are assumed to be equal
#' \item \code{theta4} is larger than \code{theta5}
#' }
#' 
#' @family functions to evaluate inequality constraints
#' @seealso \code{\link{generate_restriction_list}}
#' @examples
#' # priors
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' 
#' # informed hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4')
#' Hr            <- c('theta1', '<',  'theta2', '<', 'theta3', '<', 'theta4')
#' 
#' results_prior  <- binom_bf_inequality(Hr=Hr, a=a, b=b, 
#' factor_levels=factor_levels, prior=TRUE, seed = 2020)
#' # corresponds to
#' cbind(exp(results_prior$logml), 1/factorial(4))
#' 
#' # alternative - if you have samples and a restriction list
#' inequalities  <- generate_restriction_list(Hr=Hr, a=a,b=b,
#' factor_levels=factor_levels)$inequality_constraints
#' prior_samples <- binom_tsampling(inequalities, niter = 2e3, 
#' prior=TRUE, seed = 2020)
#' results_prior <- binom_bf_inequality(prior_samples, inequalities, seed=2020)
#' cbind(exp(results_prior$logml), 1/factorial(4))
#' 
#' @family functions to evaluate informed hypotheses
#' 
#' @references 
#' \insertRef{gronau2017tutorial}{multibridge} 
#' 
#' \insertRef{sarafoglou2020evaluatingPreprint}{multibridge} 
#' 
#' @export
binom_bf_inequality <- function(samples=NULL, restrictions=NULL,
                              x = NULL, n = NULL, Hr=NULL,
                              a = rep(1,ncol(samples)),
                              b = rep(1,ncol(samples)),
                              factor_levels = NULL,
                              prior = FALSE,
                              index = 1, maxiter = 1e3, seed=NULL,
                              niter=5e3, nburnin=niter * 0.05){
  
  # Step 1: Check for restriction list; create one if necessary
  if(is.null(restrictions) | !inherits(restrictions, 'bmult_rl') & !inherits(restrictions, 'bmult_rl_ineq')){
    
    if(is.null(factor_levels)){
      
      if(!is.null(colnames(samples))){
        
        factor_levels <- colnames(samples)
        
      } else {
        
        factor_levels <- paste0('theta', 1:ncol(samples))
        
      }
      
    }
    
    # transform 2-dimensional table to vector of counts and total
    userInput     <- .checkIfXIsVectorOrTable(x, n)
    x             <- userInput$counts
    total         <- userInput$total
    factor_levels <- .checkFactorLevels(x, factor_levels)
    .checkAlphaAndData(alpha = a, beta=b, counts = x, total=total)
    .checkNrParameters(factor_levels, alpha = a, counts = x)
    Hr            <- .checkSpecifiedConstraints(Hr, factor_levels)
    
    # Put factor levels in order for analysis
    # Put factor levels in order for analysis
    constrained_factors   <- purrr::keep(Hr, function(x) any(x %in% factor_levels))
    
    # Convert alpha vector and data vector accordingly &
    # discard data and concentration parameters from unconstrained factors
    match_sequence <- match(constrained_factors, factor_levels)
    a              <- a[match_sequence]
    b              <- b[match_sequence]
    x              <- x[match_sequence]
    total          <- total[match_sequence]
    
    restriction_list <- generate_restriction_list(Hr=Hr, factor_levels=factor_levels, a=a, b=b, x=x, n=total)
    restrictions     <- restriction_list$inequality_constraints
    
  } else {
    
    .checkRestrictionListClass(restrictions)
    
    if(inherits(restrictions, 'bmult_rl')){
      
      restrictions <- restrictions$inequality_constraints
      
    }
    
  }
  
  # Step 2: get samples if necessary
  if(!is.null(samples)){
    
    .checksIfMatrix(samples)
    
  } else {
    
    samples <- binom_tsampling(restrictions, index=index, prior=prior,
                                      seed=seed,niter=niter, nburnin=nburnin)
    
  }
  
  # Step 3: Extract Relevant Information To Start Bridge Sampling Routine
  boundaries       <- restrictions$boundaries[[index]]
  binom_equal      <- restrictions$mult_equal[[index]]
  hyp_direction    <- restrictions$direction[index]
  hyp              <- restrictions$hyp[[index]]
  
  if(prior | is.null(restrictions$counts_inequalities[[index]])){
    
    a   <- restrictions$alpha_inequalities[[index]]
    b   <- restrictions$beta_inequalities[[index]]
    
  } else if(!is.null(restrictions$counts_inequalities[[index]])) {
    
    a   <- restrictions$alpha_inequalities[[index]] + restrictions$counts_inequalities[[index]]
    b   <- restrictions$beta_inequalities[[index]] + (restrictions$total_inequalities[[index]] - restrictions$counts_inequalities[[index]])
    
  }
  
  # set seed if wanted
  if(!is.null(seed) & is.numeric(seed)){
    set.seed(seed)
  }
  
  # check if correct number of parameters were provided
  .checkNrParameters(samples = samples, boundaries = boundaries)
  
  # Step 4: Bridge Sampling
  
  ###    Code by Gronau et al. (2017) - online appendix ###
  ###    Modified by Alexandra Sarafoglou               ###
  
  # Note that before applying this function the user needs to:
  # 1. Collect 2*N1 samples from the truncated prior and posterior distribution
  #    (e.g., through MCMC sampling)
  # 2. Choose a suitable proposal distribution. Here we choose the multivariate normal &
  #    Specify the function for evaluating the log of the unnormalized density.
  #    This function is here referred to as log_unnormalized_density.
  
  # 2. Specify the function for evaluating the log of the unnormalized density
  # 3. Transform the parameters to the real line
  samples     <- tbinom_trans(samples, boundaries, binom_equal, hyp_direction)
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
  q11 <- log_unnormalized_tbinom(samples_4_iter, boundaries, binom_equal, hyp_direction, a, b)
  q21 <- log_unnormalized_tbinom(gen_samples   , boundaries, binom_equal, hyp_direction, a, b)
  
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
