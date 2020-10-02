#' @title Creates Restriction List Based On User Specified Informed Hypothesis
#'
#' @description This function encodes the user specified informed hypothesis. It creates a separate restriction list
#' for the full model, and all independent equality and inequality constraints. The returned list features relevant 
#' information for the transformation and sampling of the model parameters, such as information about the upper and 
#' lower bound for each parameter, and the indeces of equality constrained and free parameters. 
#'
#' @inheritParams binomBfInformed
#' @param a numeric. Vector with concentration parameters of Dirichlet distribution (for multinomial models) or alpha 
#' parameters for independent beta distributions (for binomial models). Default sets all parameters to 1
#' @param x a vector with data (for multinomial models) or a vector of counts of successes, or a two-dimensional table 
#' (or matrix) with 2 columns, giving the counts of successes and failures, respectively (for binomial models).
#' 
#' @return Restriction list containing the following elements:
#' \describe{
#' \item{\code{$full_model}}{\itemize{
#'   \item  \code{hyp}: character vector containing the informed hypothesis as specified by the user.
#'   \item  \code{parameters_full}: character vector containing the names for each constrained parameter.
#'   \item \code{alpha_full}: numeric vector containing the concentration parameters 
#'    of the Dirichlet distribution (when evaluating ordered multinomial parameters) or alpha parameters of the 
#'    beta distribution (when evaluating ordered binomial parameters).
#'  \item  \code{beta_full}: numeric vector containing the values of beta parameters of the beta distribution (when evaluating ordered binomial parameters).
#'  \item  \code{counts_full}: numeric vector containing data values (when evaluating multinomial parameters), or number of successes (when evaluating ordered binomial parameters).
#'  \item  \code{total_full}: numeric vector containing the number of observations (when evaluating ordered binomial 
#'  parameters, that is, number of successes and failures).
#' }}
#' \item{\code{$equality_constraints}}{\itemize{
#'   \item  \code{hyp}: list containing all independent equality constrained hypotheses.
#'   \item  \code{parameters_equality}: character vector containing the names for each equality constrained parameter.
#'   \item \code{equality_hypotheses}: list containing the indices of each equality constrained parameter. Note that these indices are based on the vector of all factor levels
#'  \item  \code{alpha_equalities}: list containing the concentration parameters for equality constrained hypotheses (when evaluating multinomial parameters) or alpha parameters of the 
#'    beta distribution (when evaluating ordered binomial parameters).
#'  \item  \code{beta_equalities}: list containing the values of beta parameters of the beta distribution (when evaluating ordered binomial parameters).
#'  \item  \code{counts_equalities}: list containing data values (when evaluating multinomial parameters), or number of successes (when evaluating ordered binomial parameters) of each equality constrained parameter.
#'  \item  \code{total_equalitiesl}: list containing the number of observations of each equality constrained parameter (when evaluating ordered binomial 
#'  parameters, that is, number of successes and failures)
#' }}
#' \item{\code{$inequality_constraints}}{\itemize{
#'   \item  \code{hyp}: list containing all independent inequality constrained hypotheses.
#'   \item  \code{parameters_inequality}: list containing the names for each inequality constrained parameter.
#'   \item \code{inequality_hypotheses}: list containing the indices of each inequality constrained parameter.
#'  \item  \code{alpha_inequalities}: list containing for inequality constrained hypotheses the concentration parameters 
#'  of the Dirichlet distribution (when evaluating ordered multinomial parameters) or alpha parameters of the beta distribution (when 
#'  evaluating ordered binomial parameters).
#'  \item  \code{beta_inequalities}: list containing for inequality constrained hypotheses the values of beta parameters of 
#'  the beta distribution (when evaluating ordered binomial parameters).
#'  \item  \code{counts_inequalities}: list containing for inequality constrained parameter data values (when evaluating 
#'  multinomial parameters), or number of successes (when evaluating ordered binomial parameters).
#'  \item  \code{total_inequalities}: list containing for each inequality constrained parameter the number of observations 
#'  (when evaluating ordered binomial parameters, that is, number of successes and failures).
#'  \item  \code{boundaries}: list that lists for each inequality constrained parameter the index of parameters that 
#'  serve as its upper and lower bounds. Note that these indices refer to the collapsed categories (i.e., categories after conditioning
#'  for equality constraints). If a lower or upper bound is missing, for instance because the current parameter is set to be the 
#'  smallest or the largest, the bounds take the value 'int(0)'.
#'  \item  \code{nr_mult_equal}: list containing multiplicative elements of collapsed categories
#'  \item  \code{nr_mult_free}: list containing multiplicative elements of free parameters.
#'  \item  \code{mult_equal}: list that contains for each lower and upper bound of each inequality constrained parameter 
#'  necessary multiplicative elements to recreate the implied order restriction, even for collapsed parameter values. If 
#'  there is no upper or lower bound, the multiplicative element will be 0.
#'  \item  \code{nineq_per_hyp}: numeric vector containing the total number of inequality constrained parameters 
#'  for each independent inequality constrained hypotheses.
#'  \item  \code{direction}: character vector containing the direction for each independent inequality constrained 
#'  hypothesis. Takes the values 'smaller' or 'larger'.
#' }}
#' }
#' @details The restriction list can be created for both binomial and multinomial models. If multinomial models are specified,
#' the arguments \code{b} and \code{n} should be left empty and \code{x} should not be a table or matrix.
#' @examples
#' # Restriction list for ordered multinomial
#' x <- c(1, 1, 1, 1)
#' a <- c(1, 1, 1, 1)
#' factor_levels <- c('mult1', 'mult2', 'mult3', 'mult4')
#' Hr <- c('mult1 > mult2 , mult3 = mult4')
#' restrictions <- generateRestrictionList(x=x, Hr=Hr, a=a, 
#' factor_levels=factor_levels)
#' @export
generateRestrictionList <- function(x=NULL, n = NULL, Hr, a, b = NULL, factor_levels) {
  
  ## Initialize Output List
  out <- list(full_model             = list(hyp             = NULL,
                                            parameters_full = NULL,
                                            alpha_full      = NULL,
                                            beta_full       = NULL,
                                            counts_full     = NULL,
                                            total_full      = NULL),
              equality_constraints   = list(hyp                 = NULL,
                                            parameters_equality = NULL,
                                            equality_hypotheses = NULL,
                                            alpha_equalities    = NULL,
                                            counts_equalities   = NULL),
              inequality_constraints = list(hyp                   = NULL,
                                            parameters_inequality = NULL,
                                            inequality_hypotheses = NULL,
                                            alpha_inequalities    = NULL,
                                            beta_inequalities     = NULL,
                                            counts_inequalities   = NULL,
                                            total_inequalities    = NULL,
                                            nr_mult_equal         = NULL,
                                            nr_mult_free          = NULL,
                                            mult_equal            = NULL,
                                            boundaries            = NULL,
                                            nineq_per_hyp         = NULL,
                                            direction             = NULL))
  
  signs <- c(equal='=', smaller='<', larger='>', free=',', linebreak='&')
  
  # transform 2-dimensional table to vector of counts and total
  userInput     <- .checkIfXIsVectorOrTable(x, n)
  x             <- userInput$counts
  n             <- userInput$total
  factor_levels <- .checkFactorLevels(x, factor_levels)
  
  # check user input
  if(is.factor(factor_levels)) factor_levels <- levels(factor_levels)
  Hr <- .checkSpecifiedConstraints(Hr, factor_levels, signs)
  
  factors_analysis <- factor_levels[factor_levels %in% Hr]
  
  ## Encode the full model
  out$full_model$hyp             <- Hr
  out$full_model$parameters_full <- factors_analysis
  out$full_model$alpha_full      <- a[which(factors_analysis %in% factor_levels)]
  out$full_model$beta_full       <- b[which(factors_analysis %in% factor_levels)]
  out$full_model$counts_full     <- x[which(factors_analysis %in% factor_levels)]
  out$full_model$total_full      <- n[which(factors_analysis %in% factor_levels)]
  
  ## Encode equality constraints: only relevant for ordered multinomial parameters
  # 2.1 saves the equality constrained hypotheses
  distinct_equalities  <- .splitAt(Hr, signs[c('free','linebreak', 'smaller', 'larger')]) %>%
    purrr::keep(function(x) any(x == signs['equal']))
  out$equality_constraints$hyp <- distinct_equalities
  # splits of categories that are not equality constrained
  equality_list <- distinct_equalities %>%
    purrr::map(function(x) which(factor_levels %in% x)) %>% purrr::compact()
  # 2.2 saves names for each factor level that is equality constrained
  out$equality_constraints$parameters_equality <- factor_levels[unlist(equality_list)]
  # 2.3 saves indices for equality constrained parameters
  out$equality_constraints$equality_hypotheses <- equality_list
  # 2.4 saves the concentration parameters for equality constrained hypotheses
  out$equality_constraints$alpha_equalities    <- a[unlist(out$equality_constraints$equality_hypotheses)]
  # 2.5 saves the number of observations per category for equality constrained hypotheses
  out$equality_constraints$counts_equalities   <- x[unlist(out$equality_constraints$equality_hypotheses)]
  
  ## Encode inequality constraints
  # 3.0 saves the free parameters
  distinct_free_parameters   <- .splitAt(Hr, signs[c('linebreak', 'smaller', 'larger')]) %>% purrr::keep(function(x) any(x == signs['free']))
  # 3.1 saves the inequality constrained hypotheses
  distinct_inequalities          <- .splitAt(Hr, signs['linebreak']) %>% purrr::keep(function(x) any(signs[c('smaller', 'larger')] %in% x))
  out$inequality_constraints$hyp <- distinct_inequalities
  
  if(!purrr::is_empty(out$inequality_constraints$hyp)){
    
    # loop over each independent order restricted hypothesis
    for(i in 1:length(distinct_inequalities)){
      
      # splits of categories that are not inequality constrained
      ineq_param_index <- which(factors_analysis %in% distinct_inequalities[[i]])
      # collapses categories that are equality constrained
      inequality_constrained_parameters <- factors_analysis[ineq_param_index]
      equalities_to_be_collapsed        <- distinct_equalities %>%
        purrr::map(function(x) which(inequality_constrained_parameters %in% x)) %>% purrr::compact()
      # 3.2 saves names for each factor level that is inequality constrained
      inequality_constrained_parameters                     <- .collapseCategories(inequality_constrained_parameters, equalities_to_be_collapsed, is_numeric_value = FALSE)
      out$inequality_constraints$parameters_inequality[[i]] <- inequality_constrained_parameters
      # 3.3 saves indices for collapsed inequalities
      out$inequality_constraints$inequality_hypotheses[[i]] <- unlist(distinct_inequalities[i] %>%
                                                                        purrr::map(function(x) which(inequality_constrained_parameters %in% x)) %>% purrr::compact())
      # 3.4 saves the concentration parameters for inequality constrained hypotheses
      alpha_inequalities                                 <- out$full_model$alpha_full[ineq_param_index]
      out$inequality_constraints$alpha_inequalities[[i]] <- .collapseCategories(alpha_inequalities, equalities_to_be_collapsed, correct = TRUE)
      
      if(!purrr::is_empty(b)){
        
        beta_inequalities                                 <- out$full_model$beta_full[ineq_param_index]
        out$inequality_constraints$beta_inequalities[[i]] <-.collapseCategories(beta_inequalities, equalities_to_be_collapsed, correct = TRUE)
        
      }
      
      if(!purrr::is_empty(x)){
        
        counts_inequalities                                 <- out$full_model$counts_full[ineq_param_index]
        out$inequality_constraints$counts_inequalities[[i]] <- .collapseCategories(counts_inequalities, equalities_to_be_collapsed)
        
      }
      
      if(!purrr::is_empty(n)){
        
        total_inequalities                                   <- out$full_model$total_full[ineq_param_index]
        out$inequality_constraints$total_inequalities[[i]]   <- .collapseCategories(total_inequalities, equalities_to_be_collapsed)
        
      }
      
      # 3.6 saves the upper and lower boundaries for each inequality constrained parameter
      # encode the direction of the hypotheses
      encode_order_direction <- distinct_inequalities[i] %>%
        purrr::map(function(x) ifelse(signs['smaller'] %in% x, 'smaller', 'larger')) %>%
        purrr::flatten_chr()
      distinct_inequality <- .splitAt(distinct_inequalities[[i]],signs[c('smaller', 'larger')])
      distinct_inequality <- rapply(distinct_inequality, function(x) which(inequality_constrained_parameters %in% x), how = 'list')
      # create numeric representation of inequality constraints
      num_representation <- sapply(distinct_inequality, length)
      
      if(encode_order_direction == 'smaller'){
        
        num_representation <- rep(1:length(num_representation), num_representation)
        
      } else {
        
        num_representation <- rep(length(num_representation):1, num_representation)
        
      }
      
      # create boundaries list for each inequality constrained parameter
      out$inequality_constraints$boundaries[[i]] <- vector('list', length(inequality_constrained_parameters))
      out$inequality_constraints$boundaries[[i]] <- lapply(out$inequality_constraints$boundaries[[i]], function(x) x <- list(lower = NULL, upper = NULL))
      
      for(j in seq_along(inequality_constrained_parameters)){
        
        current_parameter                                     <- num_representation[j]
        out$inequality_constraints$boundaries[[i]][[j]]$lower <- which(num_representation < current_parameter)
        out$inequality_constraints$boundaries[[i]][[j]]$upper <- which(num_representation > current_parameter)
      }
      
      # 3.7 create list for multiplicative elements of collapsed or free parameters
      multiplicative_elements                       <- .collapseCategories(rep(1, length(factors_analysis[ineq_param_index])), equalities_to_be_collapsed)
      out$inequality_constraints$nr_mult_equal[[i]] <- multiplicative_elements
      
      free_elements     <- distinct_free_parameters %>% purrr::map(function(x) which(inequality_constrained_parameters %in% x)) %>% purrr::compact()
      nr_mult_free     <- rep(0, length(inequality_constrained_parameters))
      
      for(j in seq_along(free_elements)) {
        
        if(encode_order_direction == 'smaller'){
          
          nr_mult_free_j                   <- sum(multiplicative_elements[free_elements[[j]]]) -
            cumsum(multiplicative_elements[free_elements[[j]]])
          nr_mult_free[free_elements[[j]]] <- nr_mult_free_j
          
        } else {
          
          nr_mult_free_j <- sum(multiplicative_elements[free_elements[[j]]]) -
            cumsum(rev(multiplicative_elements[free_elements[[j]]]))
          nr_mult_free[free_elements[[j]]] <- rev(nr_mult_free_j)
          
        }
      }
      
      out$inequality_constraints$nr_mult_free[[i]]     <- nr_mult_free
      
      # 3.8 create list of the multiplicative elements for each lower and upper bound of each inequality constrained parameter
      out$inequality_constraints$mult_equal[[i]] <- vector('list', length(inequality_constrained_parameters))
      out$inequality_constraints$mult_equal[[i]] <- lapply(out$inequality_constraints$boundaries[[i]], function(x) x <- list(lower = 1, upper = 1))
      
      for(j in seq_along(inequality_constrained_parameters)){
        
        current_mult_index   <- multiplicative_elements[j]
        lower_parameters     <- out$inequality_constraints$boundaries[[i]][[j]]$lower
        upper_parameters     <- out$inequality_constraints$boundaries[[i]][[j]]$upper
        out$inequality_constraints$mult_equal[[i]][[j]]$lower <- 1/multiplicative_elements[lower_parameters] * current_mult_index
        out$inequality_constraints$mult_equal[[i]][[j]]$upper <- 1/multiplicative_elements[upper_parameters] * current_mult_index
      }
      
      # 3.9 saves the total number of inequality constrained parameters for each independent constraint
      out$inequality_constraints$nineq_per_hyp[i] <- length(out$inequality_constraints$inequality_hypotheses[[i]])
      # 3.10 saves the direction for each independent constraint
      out$inequality_constraints$direction[i] <- encode_order_direction
    }
  }
  
  # assign class
  class(out)                        <-  append(class(out)                       , 'bmult_rl')
  class(out$full_model)             <-  append(class(out$full_model)            , 'bmult_rl_full')
  class(out$equality_constraints)   <-  append(class(out$equality_constraints)  , 'bmult_rl_eq')
  class(out$inequality_constraints) <-  append(class(out$inequality_constraints), 'bmult_rl_ineq')
  # return output
  return(out)
}
