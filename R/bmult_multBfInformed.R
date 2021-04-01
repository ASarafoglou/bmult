#' @useDynLib multibridge
NULL

#' @importFrom Rcpp evalCpp
NULL

#' @title Evaluates Informed Hypotheses on Multinomial Parameters
#' 
#' @description Evaluates informed hypotheses on multinomial parameters. These hypotheses can contain
#' (a mixture of) inequality constraints, equality constraints, and free parameters.
#' Informed hypothesis \eqn{H_r} states that category proportions obey the particular constraint.
#' \eqn{H_r} can be tested against the encompassing hypothesis \eqn{H_e} or the null hypothesis \eqn{H_0}. 
#' Encompassing hypothesis \eqn{H_e} states that category proportions are free to vary.
#' Null hypothesis \eqn{H_0} states that category proportions are exactly equal.
#' 
#'
#' @inherit mult_tsampling
#' @inherit mult_bf_inequality
#' @param x numeric. Vector with data
#' @param Hr string or character. Encodes the user specified informed hypothesis. Use either specified \code{factor_levels}
#' or indices to refer to parameters. See ``Note'' section for details on how to formulate informed hypotheses 
#' @param factor_levels character. Vector with category names. Must be the same length as \code{x}
#' @param a numeric. Vector with concentration parameters of Dirichlet distribution. Must be the same length as \code{x}. Default sets all concentration parameters to 1
#' @param cred_level numeric. Credible interval for the posterior point estimates. Must be a single number between 0 and 1
#' @param niter numeric. Vector with number of samples to be drawn from truncated distribution
#' @param bf_type character. The Bayes factor type. When the informed hypothesis is compared to the encompassing hypothesis, 
#' the Bayes factor type can be \code{LogBFer}, \code{BFer}, or \code{BFre}. When the informed hypothesis is compared to the null hypothesis, 
#' the Bayes factor type can be \code{LogBFr0}, \code{BF0r}, or \code{BFr0}. Default is \code{LogBFer}
#' @param seed numeric. Sets the seed for reproducible pseudo-random number generation
#' 
#' @details  
#' The model assumes that data follow a multinomial distribution and assigns a Dirichlet distribution as prior for the model parameters 
#' (i.e., underlying category proportions). That is:
#' \deqn{x ~ Multinomial(N, \theta)}
#' \deqn{\theta ~ Dirichlet(\alpha)}
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
#' @return List consisting of the following elements 
#' \describe{
#' \item{\code{$bf_list}}{gives an overview of the Bayes factor analysis:
#' \itemize{
#' \item \code{bf_type}: string. Contains Bayes factor type as specified by the user
#' \item \code{bf}: data.frame. Contains Bayes factors for all Bayes factor types
#' \item \code{error_measures}: data.frame. Contains for the overall Bayes factor
#' the approximate relative mean-squared error \code{re2}, the approximate coefficient of variation \code{cv}, and the approximate percentage error \code{percentage}
#' \item \code{logBFe_equalities}: data.frame. Lists the log Bayes factors for all independent equality constrained hypotheses
#' \item \code{logBFe_inequalities}: data.frame. Lists the log Bayes factor for all independent inequality constrained hypotheses
#' }}
#' \item{\code{$cred_level}}{numeric. User specified credible interval}
#' \item{\code{$restrictions}}{list that encodes informed hypothesis for each independent restriction:
#' \itemize{
#' \item \code{full_model}: list containing the hypothesis, parameter names, data and prior specifications for the full model.
#' \item \code{equality_constraints}: list containing the hypothesis, parameter names, data and prior specifications for each equality constrained hypothesis.
#' \item \code{inequality_constraints}: list containing the hypothesis, parameter names, data and prior specifications for each inequality constrained hypothesis. 
#' In addition, in \code{nr_mult_equal} and \code{nr_mult_free} encodes which and how many parameters are
#' equality constraint or free, in \code{boundaries} includes the boundaries of each parameter, in \code{nineq_per_hyp} states the number of inequality constraint 
#' parameters per independent inequality constrained hypothesis, and in \code{direction} states the direction of 
#' the inequality constraint.
#' }}
#' \item{\code{$bridge_output}}{list containing output from bridge sampling function:
#' \itemize{
#' \item \code{eval}: list containing the log prior or posterior evaluations
#' (\code{q11}) and the log proposal evaluations (\code{q12}) for the prior or posterior samples, 
#' as well as the log prior or posterior evaluations (\code{q21}) and the log proposal evaluations (\code{q22}) 
#' for the samples from the proposal distribution
#' \item \code{niter}: number of iterations of the iterative updating scheme
#' \item \code{logml}: estimate of log marginal likelihood
#' \item \code{hyp}: evaluated inequality constrained hypothesis
#' \item \code{error_measures}: list containing in \code{re2} the approximate 
#' relative mean-squared error for the marginal likelihood estimate, in \code{cv} the approximate 
#' coefficient of variation for the marginal likelihood estimate (assumes that bridge estimate is unbiased), and
#' in \code{percentage} the approximate percentage error of the marginal likelihood estimate
#' }}
#' \item{\code{$samples}}{list containing a list for prior samples and a list
#' of posterior samples from truncated distributions which were used to evaluate inequality constraints. 
#' Prior and posterior samples of independent inequality constraints are again saved 
#' in separate lists. Samples are stored as matrix of dimension \code{nsamples x nparams}.}
#' }
#' 
#' @examples 
#' # data
#' x <- c(3, 4, 10, 11, 7, 30)
#' # priors
#' a <- c(1, 1, 1, 1, 1, 1)
#' # restricted hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4', 'theta5', 
#' 'theta6')
#' Hr            <- c('theta1', '<',  'theta2', '&', 'theta3', '=', 'theta4', 
#' ',', 'theta5', '<', 'theta6')
#' output_total  <- mult_bf_informed(x, Hr, a, factor_levels, seed=2020, niter=2e3)
#' 
#' @references 
#' \insertRef{damien2001sampling}{multibridge}
#' 
#' \insertRef{gronau2017tutorial}{multibridge} 
#' 
#' \insertRef{fruhwirth2004estimating}{multibridge}
#' 
#' \insertRef{sarafoglou2020evaluatingPreprint}{multibridge} 
#' @family functions to evaluate informed hypotheses
#' @export
mult_bf_informed <- function(x, Hr, a=rep(1, length(x)), factor_levels=NULL, cred_level = 0.95, niter = 5e3, bf_type = 'LogBFer', seed=NULL, 
                           maxiter=1e3, nburnin=niter * 0.05){
  
  #######################
  ## Checks User Input ##
  #######################
  
  # transform 2-dimensional table to vector of counts and total
  bf_type   <- match.arg(bf_type, c('BFer', 'BFre', 'LogBFer', 
                                    'BF0r', 'BFr0', 'LogBFr0')) 
 
  factor_levels <- .checkFactorLevels(x, factor_levels)
  .checkCredLevel(cred_level = cred_level)
  .checkAlphaAndData(alpha = a, counts = x)
  .checkNrParameters(factor_levels, alpha = a, counts = x)
  Hr <- .checkSpecifiedConstraints(Hr, factor_levels)

  ################################
  ## Preprocessing for Analysis ##
  ################################
  
  # Put factor levels in order for analysis
  constrained_factors <- purrr::keep(Hr, function(x) any(x %in% factor_levels))
  
  # Convert alpha vector and data vector accordingly &
  # discard data and concentration parameters from unconstrained factors
  match_sequence <- match(constrained_factors, factor_levels)
  # match the order of parameters to the restricted hypothesis
  a              <- a[match_sequence]
  x              <- x[match_sequence]
  
  # Encode H_r
  restrictions          <- generate_restriction_list(Hr=Hr, factor_levels=constrained_factors, a=a, x=x)
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
      counts_equalities  <- x[equalities$equality_hypotheses[[i]]]
      thetas             <- rep(1/K_equalities, K_equalities)
      
      # conduct multinomial test for each equality constraint
      equalities_list[[i]] <- mult_bf_equality(x=counts_equalities, a=alphas_equalities, p=thetas)
      logBFe_equalities[i] <- equalities_list[[i]]$bf[['LogBFe0']]

    }
    
    logBFe_equalities <- as.data.frame(logBFe_equalities)
    
  }

  ### Evaluate inequality constraints ###
  logBFe_inequalities  <- logml_prior <- logml_post <- 0
  bs_results           <- vector('list', length(inequalities$inequality_hypotheses))
  error_measures_prior <- error_measures_post <- 0

  if(!purrr::is_empty(inequalities$hyp)){
    
    prior_samples <- post_samples <- vector('list', length(inequalities$inequality_hypotheses))
    error_measures_prior <- error_measures_post <- rep(0, length(inequalities$inequality_hypotheses))
    
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
        prior_samples[[i]]           <- mult_tsampling(inequalities, index, niter, prior=TRUE, seed=seed, nburnin=nburnin)
        colnames(prior_samples[[i]]) <- colnames_samples
        bs_results[[i]]$prior        <- mult_bf_inequality(prior_samples[[i]], restrictions=inequalities, index=index, prior=TRUE, seed=seed, maxiter=maxiter)
        logml_prior[i]               <- bs_results[[i]]$prior$logml
        error_measures_prior[i]      <- bs_results[[i]]$prior$error_measures$re2
      }

      # posterior
      post_samples[[i]]           <- mult_tsampling(inequalities, index, niter, seed=seed, nburnin=nburnin)
      colnames(post_samples[[i]]) <- colnames_samples
      bs_results[[i]]$post        <- mult_bf_inequality(post_samples[[i]], restrictions=inequalities, index=index, seed=seed, maxiter=maxiter)
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
    
    bf0_table <-  mult_bf_equality(x, a)$bf
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
  output <- list(bf_list         = bf_list,
                 cred_level      = cred_level,
                 restrictions    = restrictions,
                 bridge_output   = bs_results
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
