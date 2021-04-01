#' @title Samples From Truncated Dirichlet Density
#' 
#' @description Based on specified inequality constraints, samples from truncated 
#' prior or posterior Dirichlet density.
#' 
#' @inherit mult_bf_informed
#' @param inequalities list that contains inequality constraints for each independent inequality constrained hypotheses. The list 
#' is created in the \code{\link{generate_restriction_list}} function
#' @param index numeric. If multiple independent inequality constraints are specified, this index determines for which 
#' inequality constraint samples should be drawn. Must be a single value. Default is 1
#' @param niter numeric. A single value specifying the number of samples. Default is set to \eqn{10,000}
#' @param prior logical. If \code{TRUE} ignores the data that are encoded in \code{inequalities} and thus samples from the 
#' prior distribution. Default is \code{FALSE}.
#' @param nburnin numeric. A single value specifying the number of burn-in samples when drawing from the truncated distribution. 
#' Minimum number of burn-in samples is 10. Default is 5% of the number of samples. Burn-in samples are removed automatically after the sampling.
#' @return matrix of dimension \code{niter * nsamples} containing prior or posterior samples from truncated Dirichlet distribution.
#' @note When equality constraints are specified in the restricted hypothesis, this function samples from the conditional 
#' Dirichlet distribution given that the equality constraints hold. 
#' 
#' Only inequality constrained parameters are sampled. Free parameters or parameters that are 
#' exclusively equality constrained will be ignored. 
#' @seealso \code{\link{generate_restriction_list}}
#' @family function to sample from truncated densities
#' @examples 
#' x <- c(200, 130, 40, 10)
#' a <- c(1, 1, 1, 1)
#' factor_levels <- c('mult1', 'mult2', 'mult3', 'mult4')
#' Hr <- c('mult1 > mult2 > mult3 > mult4')
#' 
#' # generate restriction list
#' inequalities <- generate_restriction_list(x=x, Hr=Hr, a=a, 
#' factor_levels=factor_levels)$inequality_constraints
#' 
#' # sample from prior distribution
#' prior_samples <- mult_tsampling(inequalities, niter = 500, prior=TRUE)
#' # sample from posterior distribution
#' post_samples <- mult_tsampling(inequalities, niter = 500)
#' @references 
#' \insertRef{damien2001sampling}{multibridge}
#'  
#' \insertRef{sarafoglou2020evaluatingPreprint}{multibridge} 
#' @export
mult_tsampling  <- function(inequalities  , index=1, niter = 1e4, prior=FALSE, nburnin = niter*.05, seed=NULL) {
  
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
  
  # extract relevant information
  if(prior){
    prior_and_data  <- inequalities$alpha_inequalities[[index]]
  } else {
    prior_and_data  <- inequalities$alpha_inequalities[[index]] + inequalities$counts_inequalities[[index]]
  }
  
  boundaries      <- inequalities$boundaries[[index]]
  nr_mult_equal   <- inequalities$nr_mult_equal[[index]]
  mult_equal      <- inequalities$mult_equal[[index]]
  
  # logical evaluations
  bounds_per_restriction       <- sapply(boundaries, function(x) is.null(x))
  lower_bounds_per_restriction <- sapply(boundaries, function(x) length(x$lower) != 0)
  
  # define 5% of samples as burn-in; minimum number of burn-in samples is 10
  nburnin  <- max(c(10, nburnin))
  post_samples <- matrix(ncol=length(prior_and_data), nrow = (niter + nburnin))
  
  # starting values of Gibbs Sampler
  K  <- length(prior_and_data)
  z  <- rgamma(K, prior_and_data, 1)
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
      
      # check for bounds
      there_are_no_bounds <- bounds_per_restriction[k]
      
      if(there_are_no_bounds){
        
        #    sample from unrestricted gamma distribution
        z[k] <- rgamma(1, prior_and_data[k], 1)
        
      } else {
        
        #    if there are bounds
        #    sample from truncated gamma distribution
        
        # initialize lower bound
        Lo <- 0
        
        # check for lower bound
        there_is_a_lower_bound <- lower_bounds_per_restriction[k]
        if(there_is_a_lower_bound){
          
          smaller_value <- boundaries[[k]]$lower
          Lo            <- max(z[smaller_value] * mult_equal[[k]]$lower)
          
        }
        
        # check for upper bound
        there_is_a_upper_bound <- length(boundaries[[k]]$upper) != 0
        upper_bound <- ifelse(there_is_a_upper_bound, {
          
          upper_value <- boundaries[[k]]$upper
          min(z[upper_value] * mult_equal[[k]]$upper)
          
        }, 0)
        
        z[k] <- truncatedSamplingSubiteration(runif(1), runif(1), -z[k], Lo, prior_and_data[k],
                                              there_is_a_upper_bound,
                                              upper_bound)
      }
      
    }
    # 4. transform Gamma to Dirichlet samples
    post_samples[iter,] <- as.numeric(z/sum(z))
    
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