#' Create Restriction List Based On The User Specified Order Restriction
#'
#' This function encodes the user specified order restriction. It creates a separate restriction list
#' for the full model, all equality constraints, and all inequality constraints. Within each restriction list, all 
#' independent restrictions are summarized. Additionally, the function encodes relevant information for the transformation 
#' and sampling of the parameters. This information concerns the upper and lower bound for each parameter and features information
#' such as the number of equality constrained parameters within each element, or whether parameters are free to vary. 
#'
#' @param OR character vector containing the user specified order restriction. Elements within the order restriction must
#' match the factor levels of the multinomial variable. Independent restrictions are be indicated with "\n".
#' @param factor_levels character vector containing the names from all parameters.
#' @param a numeric vector with values of concentration parameters of the Dirichlet distribution (when evaluating ordered multinomial parameters) 
#' or alpha parameters of the beta distribution (when evaluating ordered binomial parameters).
#' @param b numeric vector with values of beta parameters of the beta distribution (when evaluating ordered binomial parameters).
#' @param counts (optional) numeric vector with data values (when evaluating multinomial parameters), or number of successes (when evaluating ordered 
#' binomial parameters)
#' @param total numeric vector total number of observations (when evaluating ordered binomial parameters, that is, number of successes and failures).
#' @param signs named vector indicating which symbol matches to the order restriction 'free', 'larger', 'smaller', 'equal',
#' and 'linebreak. Default is '==', '<', '>', ',', and '\n', respectively.
#' @return restriction list containing the following elements:
#'
#' (1) $full_model:
#'  1.1 sublist 'hyp': character vector containing the full ordinal restriction as specified by the user.
#'  1.2 sublist 'parameters_full': character vector containing the names for each constrained parameter.
#'  1.3 sublist 'alpha_full': numeric vector containing the concentration parameters 
#'  of the Dirichlet distribution (when evaluating ordered multinomial parameters) or alpha parameters of the 
#'  beta distribution (when evaluating ordered binomial parameters).
#'  1.4 sublist 'beta_full': numeric vector containing the values of beta parameters of 
#'  the beta distribution (when evaluating ordered binomial parameters).
#'  1.5 sublist 'counts_full': numeric vector containing data values (when evaluating 
#'  multinomial parameters), or number of successes (when evaluating ordered binomial parameters).
#'  1.6 sublist 'total_full': numeric vector containing the number of observations (when evaluating ordered binomial 
#'  parameters, that is, number of successes and failures).
#'
#' (2) $equality_constraints:
#'  2.1 sublist 'hyp': character vector containing the all independent equality constrained hypotheses.
#'  2.2 sublist 'parameters_equality': character vector containing the names for each equality constrained parameter.
#'  2.3 sublist 'equality_hypotheses': numeric vector containing the indices of each equality constrained parameter.
#'  Note that these indices are based on the vector of all factor levels
#'  2.4 sublist 'alpha_equalities': numeric vector containing the concentration parameters for equality constrained hypotheses.
#'  2.5 sublist 'counts_equalities': numeric vector containing data values for each for equality constrained parameter.
#'
#' (3) $inequality_constraints:
#'  3.1 sublist 'hyp': character vector containing the all independent inequality constrained hypotheses.
#'  3.2 sublist 'parameters_inequality': character vector containing the names for each inequality constrained parameter.
#'  3.3 sublist 'inequality_hypotheses': numeric vector containing the indices of each equality constrained parameter.
#'  3.4 sublist 'alpha_inequalities': numeric vector containing for inequality constrained hypotheses the concentration parameters 
#'  of the Dirichlet distribution (when evaluating ordered multinomial parameters) or alpha parameters of the beta distribution (when 
#'  evaluating ordered binomial parameters).
#'  3.5 sublist 'beta_inequalities': numeric vector containing for inequality constrained hypotheses the values of beta parameters of 
#'  the beta distribution (when evaluating ordered binomial parameters).
#'  3.6 sublist 'counts_inequalities': numeric vector containing for inequality constrained parameter data values (when evaluating 
#'  multinomial parameters), or number of successes (when evaluating ordered binomial parameters).
#'  3.7 sublist 'total_inequalities': numeric vector containing for each inequality constrained parameter the number of observations 
#'  (when evaluating ordered binomial parameters, that is, number of successes and failures).
#'  3.8 sublist 'boundaries': list that contains for each inequality constrained parameter the index from parameters that 
#'  serve as upper and lower bounds. Note that these indices refer to the collapsed categories. If a lower or upper bound is missing, 
#'  for instance because the current parameter is the smallest or the largest, the bounds take the value 'int(0)'.
#'  3.9 sublist 'nr_mult_equal': numeric vector containing multiplicative elements of collapsed parameters.
#'              'nr_mult_free': numeric vector containing multiplicative elements of free parameters.
#'  3.10 sublist 'mult_equal_adj': list that contains for each lower and upper bound of each inequality constrained parameter 
#'  necessary multiplicative elements to recreate the implied order restriction, even for collapsed parameter values. If 
#'  there is no upper or lower bound, the multiplicative element will be 0.
#'  3.11  sublist 'nineq_per_hyp': numeric vector containing the total number of inequality constrained parameters 
#'  for each independent inequality constrained hypotheses.
#'  3.12 sublist 'direction': character vector containing the direction for each independent inequality constrained 
#'  hypothesis. Takes the values 'smaller' or 'larger'.
#'  @export
generateRestrictionList <- function(OR, factor_levels, a, b = NULL, counts=NULL, total = NULL) {
  
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
  
  # check user input
  if(is.factor(factor_levels)) factor_levels <- levels(factor_levels)
  OR <- .checkSpecifiedConstraints(OR, factor_levels, signs)
  
  factors_analysis <- factor_levels[factor_levels %in% OR]
  
  ## Encode the full model
  out$full_model$hyp             <- OR
  out$full_model$parameters_full <- factors_analysis
  out$full_model$alpha_full      <- a[which(factors_analysis %in% factor_levels)]
  out$full_model$beta_full       <- b[which(factors_analysis %in% factor_levels)]
  out$full_model$counts_full     <- counts[which(factors_analysis %in% factor_levels)]
  out$full_model$total_full      <- total[which(factors_analysis %in% factor_levels)]
  
  ## Encode equality constraints: only relevant for ordered multinomial parameters
  # 2.1 saves the equality constrained hypotheses
  distinct_equalities  <- .splitAt(OR, signs[c('free','linebreak', 'smaller', 'larger')]) %>%
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
  out$equality_constraints$counts_equalities   <- counts[unlist(out$equality_constraints$equality_hypotheses)]
  
  ## Encode inequality constraints
  # 3.0 saves the free parameters
  distinct_free_parameters   <- .splitAt(OR, signs[c('linebreak', 'smaller', 'larger')]) %>% purrr::keep(function(x) any(x == signs['free']))
  # 3.1 saves the inequality constrained hypotheses
  distinct_inequalities          <- .splitAt(OR, signs['linebreak']) %>% purrr::keep(function(x) any(signs[c('smaller', 'larger')] %in% x))
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
      alpha_inequalities            <- out$full_model$alpha_full[ineq_param_index]
      # only relevant for ordered multinomial parameters: collpased categories
      alpha_inequalities_collapsed  <- .collapseCategories(alpha_inequalities, equalities_to_be_collapsed, correct = TRUE)
      out$inequality_constraints$alpha_inequalities[[i]] <- alpha_inequalities_collapsed
      
      if(!purrr::is_empty(b)){
        
        out$inequality_constraints$beta_inequalities[[i]] <- out$full_model$beta_full[ineq_param_index]
        
      }
      
      if(!purrr::is_empty(counts)){
        
        # 3.5 saves the number of observations for inequality constrained hypotheses
        counts_inequalities                                 <- out$full_model$counts_full[ineq_param_index]
        out$inequality_constraints$counts_inequalities[[i]] <- .collapseCategories(counts_inequalities, equalities_to_be_collapsed)
        
      }
      
      if(!purrr::is_empty(total)){
        
        out$inequality_constraints$total_inequalities[[i]]   <- out$full_model$total_full[ineq_param_index]
        
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
