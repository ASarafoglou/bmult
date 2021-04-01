#' @title Samples From Truncated Beta Densities
#' 
#' @description Based on specified inequality constraints, samples from truncated 
#' prior or posterior beta densities.
#'
#' @inherit binom_bf_informed
#' @inheritParams binom_bf_informed
#' @inheritParams mult_tsampling
#' @note When equality constraints are specified in the restricted hypothesis, this function samples from the conditional 
#' Beta distributions given that the equality constraints hold. 
#' 
#' Only inequality constrained parameters are sampled. Free parameters or parameters that are 
#' exclusively equality constrained will be ignored. 
#' @return matrix of dimension \code{niter * nsamples} containing samples from truncated beta distributions.
#' @examples 
#' x <- c(200, 130, 40, 10)
#' n <- c(200, 200, 200, 200)
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' factor_levels <- c('binom1', 'binom2', 'binom3', 'binom4')
#' Hr <- c('binom1 > binom2 > binom3 > binom4')
#' 
#' # generate restriction list
#' inequalities <- generate_restriction_list(x=x, n=n, Hr=Hr, a=a, b=b, 
#' factor_levels=factor_levels)$inequality_constraints
#' 
#' # sample from prior distribution
#' prior_samples <- binom_tsampling(inequalities, niter = 500, 
#' prior=TRUE)
#' # sample from posterior distribution
#' post_samples <- binom_tsampling(inequalities, niter = 500)
#' @seealso \code{\link{generate_restriction_list}}
#' @family functions to sample from truncated densities
#' @references 
#' \insertRef{damien2001sampling}{multibridge} 
#' @export
binom_tsampling  <- function(inequalities, index=1, niter = 1e4, prior=FALSE, nburnin = niter*.05, seed=NULL) {
  
  # if order restriction is given as character vector, create restriction list
  if(inherits(inequalities, 'bmult_rl') | inherits(inequalities, 'bmult_rl_ineq')){
    
    if(inherits(inequalities, 'bmult_rl')){
      
      inequalities  <- inequalities$inequality_constraints
      
    } 
    
  } else {
    
    stop('Provide a valid restriction list. The restriction list needs to be an object of class bmult_rl or bmult_rl_ineq as returned from generate_restriction_list.')
    
  }
  
  # set precision for large numbers
  initBigNumbers(200)
  # set seed if wanted
  if(!is.null(seed) & is.numeric(seed)){
    set.seed(seed)
  }
  
  # make sure we are in the binomial case
  .checksIfBinomial(beta=inequalities$beta_inequalities[[index]], 
                    counts=inequalities$counts_inequalities[[index]], 
                    total=inequalities$total_inequalities[[index]])
  
  # extract relevant information
  if(prior){
    a  <- inequalities$alpha_inequalities[[index]]
    b  <- inequalities$beta_inequalities[[index]]
  } else {
    a <- inequalities$alpha_inequalities[[index]] + inequalities$counts_inequalities[[index]]
    b <- inequalities$beta_inequalities[[index]] + (inequalities$total_inequalities[[index]] - inequalities$counts_inequalities[[index]])
  }
  
  boundaries      <- inequalities$boundaries[[index]]
  
  # logical evaluations
  bounds_per_restriction       <- sapply(boundaries, function(x) is.null(x))
  lower_bounds_per_restriction <- sapply(boundaries, function(x) length(x$lower) != 0)
  upper_bounds_per_restriction <- sapply(boundaries, function(x) length(x$upper) != 0)
  
  # define 5% of samples as burn-in; minimum number of burn-in samples is 10
  nburnin      <- max(c(10, nburnin))
  K            <- length(a)
  post_samples <- matrix(ncol=K, nrow = (niter + nburnin))
  
  # starting values of Gibbs Sampler
  theta <- stats::rbeta(K, a, b)
  
  iteration <- 0
  
  # initialize progress bar
  if(prior){
    progress_bar_text  <- paste0("restr. ", index , ".     prior sampling completed: [:bar] :percent time remaining: :eta")
  } else {
    progress_bar_text  <- paste0("restr. ", index , ". posterior sampling completed: [:bar] :percent time remaining: :eta")
  }
  pb <- progress::progress_bar$new(format = progress_bar_text, total = 100, clear = FALSE, width= 80)
  
  for(iter in 1:(niter+nburnin)){
    
    for(k in 1:K){
      
      smaller_value          <- boundaries[[k]]$lower
      larger_value           <- boundaries[[k]]$upper
      # check for bounds
      there_are_no_bounds    <- bounds_per_restriction[k]
      there_is_a_lower_bound <- lower_bounds_per_restriction[k]
      there_is_a_upper_bound <- upper_bounds_per_restriction[k]
      
      if(there_are_no_bounds){
        
        #    sample from unrestricted gamma distribution
        theta[k] <- stats::rbeta(1, a[k], b[k])
        
      } else {
        # if there are bounds sample from truncated beta distribution
        
        # if beta == 1; no need to introduce latent variable y
        
        if(b[k] == 1){
          
          # upper and lower bound for x
          l_theta     <- ifelse(there_is_a_lower_bound, max(theta[smaller_value]), 0)
          u_theta     <- ifelse(there_is_a_upper_bound, min(theta[larger_value]) , 1)
          
          # inverse CDF technique: (((runif(1) * (u_theta^alpha - l_theta^alpha)) + l_theta^alpha))^(1/alpha)
          theta[k] <- truncatedSamplingSubiterationBinomialCDF(runif(1), a[k], l_theta, u_theta)
          
        } else {
          
          # latent variable y: runif(1) * ((1 - theta[k])^(beta - 1))
          
          beta_minus_one <- b[k] - 1
          theta_bound    <- truncatedSamplingSubiterationBinomialY(runif(1), theta[k], beta_minus_one)
          
          # upper and lower bound for y
          l_theta     <- ifelse(there_is_a_lower_bound, max(theta[smaller_value]), 0)
          u_theta     <- ifelse(there_is_a_upper_bound, min(theta[larger_value]), 1)
          
          l_y          <- ifelse(b[k] > 1, l_theta, max(l_theta, theta_bound))
          u_y          <- ifelse(b[k] > 1, min(u_theta, theta_bound), u_theta)
          
          # inverse CDF technique: (((runif(1) * (u_y^alpha - l_y^alpha)) + l_y^alpha))^(1/alpha)
          theta[k] <- truncatedSamplingSubiterationBinomialCDF(runif(1), a[k], l_y, u_y)
          
        }
      }
      
    }
    # 4. transform Gamma to Dirichlet samples
    post_samples[iter,] <- theta
    
    # show progress
    if (iter %in% (niter/100 * seq(1, 100))) {
      pb$tick()
    }
    # if (iter %in% (niter/100 * seq(1, 100, by = 10))) {
    #   iteration <- iteration + 10
    #   print(paste('sampling completed:', iteration, '%', collapse = '\n'))
    # }
  }
  post_samples <- post_samples[-(1:nburnin), ]
  return(post_samples)
}
