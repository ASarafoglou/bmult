#' @useDynLib multibridge
NULL

#' @importFrom Rcpp evalCpp
NULL

#' @title Evaluates Informed Hypotheses on Multinomial Parameters
#' 
#' @description Evaluates Informed Hypotheses on Multinomial Parameters using bridge sampling.
#' Restricted hypothesis \eqn{H_r} states that category proportions obey a particular constraint.
#' Alternative hypothesis \eqn{H_e} states that category proportions are free to vary.
#'
#' @param x numeric. Vector with data
#' @param Hr character. Vector encoding the user specified restricted hypothesis. Use either specified factor_levels or indices to refer to parameters
#' @param factor_levels character. Vector with category names
#' @param a numeric. Vector with concentration parameters of Dirichlet distribution. Default sets all concentration parameters to 1
#' @param cred_level numeric. Credible interval for the posterior point estimates. Must be a single number between 0 and 1
#' @param niter numeric. Vector with number of samples to be drawn
#' @param bf_type character. The Bayes factor type. Ean be either 'LogBFer', 'BFer', or 'BFre'. Default is 'LogBFer'
#' @param seed numeric. Sets the seed for version control
#' 
#' @details The model assumes that data follow a multinomial distribution and assigns a Dirichlet distribution as prior for the model parameters 
#' (i.e., underlying category proportions). That is:
#' \deqn{x ~ Multinomial(N, \theta)}
#' \deqn{\theta ~ Dirichlet(\alpha)}
#' 
#' The following signs can be used to encode restricted hypotheses: \code{"<"} and \code{">"} for inequality constraints, \code{=} for equality constraints,
#' \code{","} for free parameters, and \code{"&"} for independent hypotheses. The restricted hypothesis can either be a string or a character vector.
#' For instance, the hypothesis \code{c("theta1 < theta2, theta3")} means 
#' \itemize{
#' \item \code{theta1} is smaller than both \code{theta2} and \code{theta3}
#' \item The parameters \code{theta2} and \code{theta3} both have \code{theta1} as lower bound, but are not influenced by each other.
#' }
#' The hypothesis \code{c("theta1 < theta2 = theta3 & theta4 > theta5")} means that 
#' \itemize{
#' \item \code{theta1} is smaller than \code{theta2} and \code{theta3}
#' \item \code{theta2} and \code{theta3} are assumed to be equal
#' \item \code{theta4} is larger than \code{theta5}
#' \item The restrictions on the parameters \code{theta1}, \code{theta2}, and \code{theta3} do
#' not influence the restrictions on the parameters \code{theta4} and \code{theta5}.
#' }
#' 
#' @return List consisting of the following elements 
#' \itemize{
#' \item \code{BF}: Bayes factor for restricted hypothesis compared to the encompassing hypothesis
#' \item \code{cred_level}: User specified credible interval
#' \item \code{restrictions}: full restriction list
#' \item \code{logBFe_equalities}: log Bayes factor for equality constrained parameters
#' \item \code{logBFe_inequalities}: log Bayes factor for inequality constrained parameters
#' }
#' 
#' @examples 
#' # data
#' x <- c(3, 4, 10, 11, 7, 30)
#' # priors
#' a <- c(1, 1, 1, 1, 1, 1)
#' # restricted hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6')
#' Hr            <- c('theta1', '<',  'theta2', '&', 'theta3', '=', 'theta4', ',', 'theta5', '<', 'theta6')
#' output_total  <- multBayesInformed(x, Hr, a, factor_levels, seed=2020)
#' 
#' @references 
#' \insertRef{sarafoglou2020evaluatingPreprint}{multibridge} 
#' @family functions to evaluate informed hypotheses
#' @export
multBayesInformed <- function(x, Hr, a=rep(1, length(x)), factor_levels=NULL, cred_level = 0.95, niter = 5e3, bf_type = 'LogBFer', seed=NULL){
  
  #######################
  ## Checks User Input ##
  #######################
  
  # transform 2-dimensional table to vector of counts and total
  counts        <- x
  
  factor_levels <- .checkFactorLevels(x, factor_levels)
  .checkCredLevel(cred_level = cred_level)
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
  counts                <- counts[match_sequence]
  
  # Encode H_r
  restrictions          <- generateRestrictionList(Hr=Hr, factor_levels=constrained_factors, a=a, x=counts)
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
  if(!purrr::is_empty(equalities$hyp)){

    for(i in seq_along(equalities$equality_hypotheses)){

      # extract relevant prior information and data
      K_equalities       <- length(equalities$equality_hypotheses[[i]])
      alphas_equalities  <- a[equalities$equality_hypotheses[[i]]]
      counts_equalities  <- counts[equalities$equality_hypotheses[[i]]]
      thetas             <- rep(1/K_equalities, K_equalities)
      
      # conduct multinomial test for each equality constraint
      equalities_list[[i]] <- multBfEquality(x=counts_equalities, a=alphas_equalities, theta=thetas)
      logBFe_equalities[i] <- equalities_list[[i]]$bf[['LogBFe0']]

    }
    
    logBFe_equalities <- as.data.frame(logBFe_equalities)
    
  }

  ### Evaluate inequality constraints ###
  logBFe_inequalities <- logml_prior <- logml_post <- 0
  bs_results <- list()

  if(!purrr::is_empty(inequalities$hyp)){

    prior.samples <- post.samples <- vector('list', length(inequalities$inequality_hypotheses))
    
    for(i in seq_along(inequalities$inequality_hypotheses)){
      
      index            <- i
      colnames_samples <- inequalities$parameters_inequality[[i]]
      
      # prior
      prior_is_uniform        <- all(inequalities$alpha_inequalities[[i]] == 1)
      any_free_parameters     <- any(stringr::str_detect(inequalities$hyp[[i]], ','))
      no_collapsed_categories <- all(inequalities$nr_mult_equal[[i]] == 1)

      if(prior_is_uniform & !any_free_parameters & no_collapsed_categories){

        K_inequalities  <- inequalities$nineq_per_hyp[i]
        logml_prior[i]  <- sum(-(lfactorial(K_inequalities)))

      } else {
        prior.samples[[i]]           <- multTruncatedSampling(inequalities, index, niter, prior=TRUE, seed=seed)
        colnames(prior.samples[[i]]) <- colnames_samples
        bs_results[[i]]              <- multBfInequality(prior.samples[[i]], restrictions=inequalities, index=index, prior=TRUE, seed=seed)
        logml_prior[i]               <- bs_results[[i]]$logml
      }

      # posterior
      post.samples[[i]]           <- multTruncatedSampling(inequalities, index, niter, seed=seed)
      colnames(post.samples[[i]]) <- colnames_samples
      bs_results[[i]]             <- multBfInequality(post.samples[[i]], restrictions=inequalities, index=index, seed=seed)
      logml_post[i]               <- bs_results[[i]]$logml
    }

    # compute BF_inequality(inequality|equality)
   logBFe_inequalities    <- logml_prior - logml_post
  }
  
  ### Compute Bayes Factor BF_er ##
  logBFer <- sum(logBFe_equalities) + sum(logBFe_inequalities)
  
  bf_list <- list(bf_type = bf_type,
                  bf      = data.frame(LogBFer = logBFer,
                                       BFer    = exp(logBFer),
                                       BFre    = 1/exp(logBFer)))
  
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
    
    output$bf_list$equalities_list   <- equalities_list
    output$bf_list$logBFe_equalities <- logBFe_equalities
    
  }
  
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
