#' Evaluates Informed Hypotheses on Multiple Binomial Parameters
#'
#' Computes Bayes factor for inequality constrained binomial parameters using the encompassing prior approach.
#' Restricted hypothesis Hr states that binomial proportions follow a particular trend.
#' Alternative hypothesis He states that binomial proportions are free to vary.
#' 
#' @param factor_levels character vector with categories
#' @param Hr character vector encoding the user specified order restriction
#' @param a numeric vector with alpha parameters
#' @param b numeric vector with beta parameters
#' @param counts numeric vector with number of successes
#' @param total numeric vector with total number of observations
#' @param niter number of samples to be drawn
#' @param bf_type 'bf_type'. The Bayes factor type can be either 'LogBFer', 'BFer', or 'BFre'. Default is 'LogBFer'.
#' @param seed set the seed for version control.
#' @param ... additional arguments (currently ignored).
#' @return list consisting of the following elements:
#'         (1) BF: Bayes factor for restricted hypothesis compared to the encompassing hypothesis
#'         (2) restrictions: full restriction list
#'         (4) logBFe_inequalities: log Bayes factor for inequality constrained parameters
#' @examples
#' # data
#' counts <- c(3, 4, 10, 11)
#' total  <- c(15, 12, 12, 12)
#' # priors
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' # restricted hypothesis
#' factor_levels <- c('binom1', 'binom2', 'binom3', 'binom4')
#' Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
#' output_total  <- binomBayesInformed(factor_levels, Hr, a, b, counts, total, niter = 5e4, bf_type = 'LogBFer', seed=2020)
#' @export
binomBayesInformed <- function(factor_levels, Hr, a, b, counts, total, niter = 5e4, bf_type = 'LogBFer', seed=NULL, ...){
  
  #######################
  ## Checks User Input ##
  #######################
  
  if(is.factor(factor_levels)) factor_levels <- levels(factor_levels)
  .checkAlphaAndData(alpha = a, counts = counts)
  .checkNrParameters(factor_levels, alpha = a, counts = counts)
  # restriction_signs <- .checkRestrictionSigns(restriction_signs)
  Hr <- .checkSpecifiedConstraints(Hr, factor_levels)
  
  ################################
  ## Preprocessing for Analysis ##
  ################################
  
  # Put factor levels in order for analysis
  constrained_factors   <- purrr::keep(factor_levels, function(x) any(x %in% Hr))
  
  # Convert alpha vector and data vector accordingly &
  # discard data and concentration parameters from unconstrained factors
  match_sequence        <- order(na.omit(match(factor_levels, constrained_factors)))
  a                     <- a[match_sequence]
  b                     <- b[match_sequence]
  counts                <- counts[match_sequence]
  total                 <- total[match_sequence]
  
  # Encode H_r
  restrictions          <- generateRestrictionList(Hr, constrained_factors, a=a, b=b, counts=counts, total=total, binom=TRUE)
  inequalities          <- restrictions$inequality_constraints
  boundaries            <- inequalities$boundaries
  ninequalities         <- inequalities$nineq_per_hyp
  
  ##############
  ## Analysis ##
  ##############
  
  ### Evaluate inequality constraints ###
  logBFe_inequalities <- logml_prior <- logml_post <- 0
  bs_results <- list()
  
  if(!purrr::is_empty(inequalities$hyp)){
    
    prior.samples <- post.samples <- vector('list', length(inequalities$inequality_hypotheses))
    
    for(i in seq_along(inequalities$inequality_hypotheses)){
      
      index            <- i
      colnames_samples <- inequalities$parameters_inequality[[i]]
      
      # prior
      prior_is_uniform        <- all(inequalities$alpha_inequalities[[i]] == 1 & inequalities$beta_inequalities[[i]] == 1)
      any_free_parameters     <- any(stringr::str_detect(inequalities$hyp[[i]], ','))
      
      if(prior_is_uniform & !any_free_parameters){
        
        K_inequalities  <- inequalities$nineq_per_hyp[i]
        logml_prior[i]  <- sum(-(lfactorial(K_inequalities)))
        
      } else {
        prior.samples[[i]]           <- binomTruncatedSampling(inequalities, index, niter, prior = TRUE, seed = seed)
        colnames(prior.samples[[i]]) <- colnames_samples
        bs_results[[i]]              <- binomBfInequality(prior.samples[[i]], restrictions=inequalities, index=index, prior=TRUE, seed=seed)
        logml_prior[i]               <- bs_results[[i]]$logml
      }
      
      # posterior
      post.samples[[i]]           <- binomTruncatedSampling(inequalities, index, niter, seed = seed)
      colnames(post.samples[[i]]) <- colnames_samples
      bs_results[[i]]             <- binomBfInequality(post.samples[[i]], restrictions=inequalities, index=index, seed=seed)
      logml_post[i]               <- bs_results[[i]]$logml
    }
    
    # compute BF_inequality(inequality|equality)
    logBFe_inequalities    <- logml_prior - logml_post
  }
  
  ### Compute Bayes Factor BF_er ##
  logBFer <- sum(logBFe_inequalities)
  # logBFer <- sum(logBFe_equalities) + sum(logBFe_inequalities)
  
  bf_list <- list(bf_type = bf_type,
                  bf      = data.frame(LogBFer = logBFer,
                                       BFer    = exp(logBFer),
                                       BFre    = 1/exp(logBFer)))
  
  ######################
  # Create Output List #
  ######################
  
  # Bayes factors
  output <- list(bf_list             = bf_list,
                 restrictions        = restrictions,
                 bridge_output       = bs_results
  )
  
  # # More information about equality constraints
  # if(!purrr::is_empty(equalities$hyp)){
  #   
  #   output$bf_list$equalities_list   <- equalities_list
  #   output$bf_list$logBFe_equalities <- logBFe_equalities
  #   
  # }
  
  # More information about inequality constraints
  if(!purrr::is_empty(inequalities$hyp)){
    
    output$samples <- list(post.samples  = post.samples,
                           prior.samples = prior.samples)
    
    output$bf_list$logBFe_inequalities  <- data.frame(
      logBFe_inequalities = logBFe_inequalities, 
      logml_prior         = logml_prior, 
      logml_post          = logml_post
    )
    
  }
  
  # assign class
  class(output) <- 'bmult'
  
  invisible(output)
}
