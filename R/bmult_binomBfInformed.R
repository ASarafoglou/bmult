#' @title Evaluates Informed Hypotheses on Multiple Binomial Parameters
#'
#' @description Evaluates informed hypotheses on multiple binomial parameters.
#' These hypotheses can contain (a mixture of) inequality constraints, equality constraints, and free parameters.
#' Informed hypothesis \eqn{H_r} states that binomial proportions obey a particular constraint.
#' \eqn{H_r} can be tested against the encompassing hypothesis \eqn{H_e} or the null hypothesis \eqn{H_0}. 
#' Encompassing hypothesis \eqn{H_e} states that binomial proportions are free to vary.
#' Null hypothesis \eqn{H_0} states that category proportions are exactly equal.
#' 
#' @inheritParams mult_bf_informed
#' @inherit mult_bf_informed 
#' @inherit binom_tsampling
#' @inherit binom_bf_inequality
#' 
#' @param x a vector of counts of successes, or a two-dimensional table (or matrix) with 2 columns, giving the counts of successes 
#' and failures, respectively
#' @param n numeric. Vector of counts of trials. Must be the same length as \code{x}. Ignored if \code{x} is a matrix or a table
#' @param a numeric. Vector with alpha parameters. Must be the same length as \code{x}. Default sets all alpha parameters to 1
#' @param b numeric. Vector with beta parameters. Must be the same length as \code{x}. Default sets all beta parameters to 1
#' 
#' @details The model assumes that the data in \code{x} (i.e., \eqn{x_1, ..., x_K}) are the observations of \eqn{K} independent
#' binomial experiments, based on \eqn{n_1, ..., n_K} observations. Hence, the underlying likelihood is the product of the 
#' \eqn{k = 1, ..., K} individual binomial functions:
#' \deqn{(x_1, ... x_K) ~ \prod Binomial(N_k, \theta_k)} 
#' Furthermore, the model assigns a beta distribution as prior to each model parameter 
#' (i.e., underlying binomial proportions). That is:
#' \deqn{\theta_k ~ Beta(\alpha_k, \beta_k)}
#' 
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
#' @examples
#' # data
#' x <- c(3, 4, 10, 11)
#' n <- c(15, 12, 12, 12)
#' # priors
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' # informed hypothesis
#' factor_levels <- c('binom1', 'binom2', 'binom3', 'binom4')
#' Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
#' output_total  <- binom_bf_informed(x, n, Hr, a, b, niter=2e3, factor_levels, seed=2020)
#' 
#' @family functions to evaluate informed hypotheses
#' @export
binom_bf_informed <- function(x, n, Hr, a, b, factor_levels=NULL, cred_level = 0.95, niter = 5e3, bf_type = 'LogBFer', seed=NULL, maxiter=1e3, nburnin=niter * 0.05){
  
  #######################
  ## Checks User Input ##
  #######################
  
  # transform 2-dimensional table to vector of counts and total
  bf_type   <- match.arg(bf_type, c('BFer', 'BFre', 'LogBFer', 
                                    'BF0r', 'BFr0', 'LogBFr0')) 
  userInput <- .checkIfXIsVectorOrTable(x, n)
  x         <- userInput$counts
  total     <- userInput$total
  
  factor_levels <- .checkFactorLevels(x, factor_levels)
  .checkCredLevel(cred_level = cred_level)
  .checkAlphaAndData(alpha = a, beta=b, counts = x, total=total)
  .checkNrParameters(factor_levels, alpha = a, counts = x)
  # restriction_signs <- .checkRestrictionSigns(restriction_signs)
  Hr <- .checkSpecifiedConstraints(Hr, factor_levels)
  
  ################################
  ## Preprocessing for Analysis ##
  ################################
  
  # Put factor levels in order for analysis
  constrained_factors   <- purrr::keep(Hr, function(x) any(x %in% factor_levels))
  
  # Convert alpha vector and data vector accordingly &
  # discard data and concentration parameters from unconstrained factors
  match_sequence <- match(constrained_factors, factor_levels)
  a              <- a[match_sequence]
  b              <- b[match_sequence]
  x              <- x[match_sequence]
  total          <- total[match_sequence]
  
  # Encode H_r
  restrictions          <- generate_restriction_list(Hr=Hr, factor_levels=constrained_factors, a=a, b=b, x=x, n=total)
  inequalities          <- restrictions$inequality_constraints
  boundaries            <- inequalities$boundaries
  ninequalities         <- inequalities$nineq_per_hyp
  equalities            <- restrictions$equality_constraints
  
  ##############
  ## Analysis ##
  ##############
  
  ### Evaluate equality constraints ###
  logBFe_equalities <- rep(0, length(equalities$hyp)) 
  equalities_list   <- list()
  # To-Do: Are Adjusted Priors Present?
  # adjustedPriorsPresent <- .checkAdjustedPriors(adjusted_priors_for_equalities, equalities$equality_hypotheses)

  if(!purrr::is_empty(equalities$hyp)){
    
    for(i in seq_along(equalities$equality_hypotheses)){
      
      # extract relevant prior information and data
      K_equalities       <- length(equalities$equality_hypotheses[[i]])
      alphas_equalities  <- a[equalities$equality_hypotheses[[i]]]
      betas_equalities   <- b[equalities$equality_hypotheses[[i]]]
      counts_equalities  <- x[equalities$equality_hypotheses[[i]]]
      total_equalities   <- total[equalities$equality_hypotheses[[i]]]
      
      # conduct multinomial test for each equality constraint
      equalities_list[[i]] <- binom_bf_equality(x=counts_equalities, n=total_equalities, a=alphas_equalities, b=betas_equalities)
      logBFe_equalities[i] <- equalities_list[[i]]$bf[['LogBFe0']]
      
    }
    
    logBFe_equalities <- as.data.frame(logBFe_equalities)
    
  }
  
  ### Evaluate inequality constraints ###
  logBFe_inequalities <- logml_prior <- logml_post <- 0
  bs_results <- vector('list', length(inequalities$inequality_hypotheses))
  error_measures_prior <- error_measures_post <- 0
  
  if(!purrr::is_empty(inequalities$hyp)){
    
    prior_samples <- post_samples <- vector('list', length(inequalities$inequality_hypotheses))
    error_measures_prior <- error_measures_post <- rep(0, length(inequalities$inequality_hypotheses))
    
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
        prior_samples[[i]]           <- binom_tsampling(inequalities, index, niter, prior = TRUE, seed = seed, nburnin=nburnin)
        colnames(prior_samples[[i]]) <- colnames_samples
        bs_results[[i]]$prior        <- binom_bf_inequality(prior_samples[[i]], restrictions=inequalities, index=index, prior=TRUE, seed=seed, maxiter=maxiter)
        logml_prior[i]               <- bs_results[[i]]$prior$logml
        error_measures_prior[i]      <- bs_results[[i]]$prior$error_measures$re2
      }
      # posterior
      post_samples[[i]]           <- binom_tsampling(inequalities, index, niter, seed = seed)
      colnames(post_samples[[i]]) <- colnames_samples
      bs_results[[i]]$post        <- binom_bf_inequality(post_samples[[i]], restrictions=inequalities, index=index, seed=seed, maxiter=maxiter)
      logml_post[i]               <- bs_results[[i]]$post$logml
      error_measures_post[i]      <- bs_results[[i]]$post$error_measures$re2
    }
    
    # compute BF_inequality(inequality|equality)
    logBFe_inequalities    <- logml_prior - logml_post
    
  }
  
  ### Compute Bayes Factor BF_er ##
  logBFer <- sum(logBFe_equalities) + sum(logBFe_inequalities)
  # compute associate error term
  re2 <- sum(error_measures_prior, error_measures_post)
  error_measures <- data.frame(re2 = re2, 
                               cv  = sqrt(re2), 
                               percentage = paste0(round(sqrt(re2)*100, 4), '%'),
                               stringsAsFactors = FALSE)
  
  if(bf_type %in% c('BF0r', 'BFr0', 'LogBFr0')){
    
    bf0_table <-  binom_bf_equality(x, n=total, a, b)$bf
    bfr_table <- data.frame(LogBFer=logBFer , 
                            BFer=exp(logBFer), 
                            BFre=1/exp(logBFer))
    
    logBFe0    <- bf0_table$LogBFe0
    logBFr0    <-  -logBFer + logBFe0
    
    bf_list <- list(bf_type    = bf_type,
                    bf         = data.frame(LogBFr0 = logBFr0,
                                            BF0r    = 1/exp(logBFr0),
                                            BFr0    = exp(logBFr0)),
                    bf0_table  = bf0_table,
                    bfr_table  = bfr_table)
    
  } else {
    
    bf_list <- list(bf_type = bf_type,
                    bf      = data.frame(LogBFer = logBFer,
                                         BFer    = exp(logBFer),
                                         BFre    = 1/exp(logBFer)))
    
  }
  
  bf_list$error_measures <- error_measures
  
  ######################
  # Create Output List #
  ######################
  
  # Bayes factors
  output <- list(bf_list             = bf_list,
                 cred_level          = cred_level,
                 restrictions        = restrictions,
                 bridge_output       = bs_results
  )
  
  # More information about equality constraints
  if(!purrr::is_empty(equalities$hyp)){
    
    output$bf_list$logBFe_equalities <- logBFe_equalities

  }
  
  # More information about inequality constraints
  if(!purrr::is_empty(inequalities$hyp)){
    
    output$samples <- list(post_samples  = post_samples,
                           prior_samples = prior_samples)
    
    output$bf_list$logBFe_inequalities  <- data.frame(
      logBFe_inequalities = logBFe_inequalities, 
      logml_prior         = logml_prior, 
      logml_post          = logml_post
    )
    
  }
  
  # assign class
  class(output) <- 'bmult'
  
  return(output)
}
