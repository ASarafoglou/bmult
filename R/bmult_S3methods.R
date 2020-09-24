#' S3 method for class 'restriction_list.bmult'
#' 
#' @export
restriction_list <- function (x, restrictions = 'inequalities') {
  UseMethod("restriction_list")
}
#' Extracts restriction list from an object of class \code{bmult}
#'
#' @param x object of class \code{bmult} as returned from \code{multBayesInformed}
#' @param restrictions specifies whether to extract restriction list for \code{equalities} or \code{inequalities. Default is \code{inequalities}.
#' @return Extracts restriction list and associated hypothesis from an object of class \code{bmult}
#' @export
restriction_list.bmult <- function(x, restrictions = 'inequalities'){
  
  if (restrictions == 'inequalities'){
    
    hyp          <- x$restrictions$inequality_constraints$hyp
    restrictions <- x$restrictions$inequality_constraints$boundaries
    expr         <- 'inequality'
    
  } else if (restrictions == 'equalities') {
    
    
    hyp          <- x$restrictions$equality_constraints$hyp
    restrictions <- x$restrictions$equality_constraints$equality_hypotheses
    expr         <- 'equality'
    
  }
  
  output <- list(hyp = hyp, restriction_list = restrictions)
  
  if(is.null(output)){
    
    cat(paste0("\nNo restriction list created for ", expr, " constrained parameters.\n\n"))
    
  } 
  
  return(output)
}

#' S3 method for class 'bridge_output.bmult'
#' 
#' @export
bridge_output <- function (x) {
UseMethod("bridge_output")
}
#' Extracts bridge sampling output from an object of class \code{bmult}
#'
#' @param x object of class \code{bmult_bridge} as returned from \code{multBayesInformed}
#' @return Extracts restriction list from an object of class \code{bmult}. The bridge sampling output contains the following elements:
#'         (1) eval: list consisting of the following elements
#'                   q11: log posterior evaluations for posterior samples.
#'                   q12: log proposal evaluations for posterior samples.
#'                   q21: log posterior evaluations for samples from proposal.
#'                   q22: log proposal evaluations for samples from proposal.
#'         (2) niter: number of iterations of the iterative updating scheme.
#'         (3) logml: estimate of log marginal likelihood
#'         (4) hyp: character vector that containts the inequality constrained hypothesis 
#'         (5) error_measures: error term of the bridge sampling estimate
#' @export
bridge_output.bmult <- function(x){
  output <- x$bridge_output
  
  if(is.null(output)){
    cat("Bridge sampling was not applied.")
  }
  
  return(output)
}

#' S3 method for class 'samples.bmult'
#' 
#' @export
samples <- function (x) {
  UseMethod("samples")
}
#' Extracts prior and posterior samples (if applicable) from an object of class \code{bmult}
#' 
#' @param x object of class \code{bmult_bridge} as returned from \code{multBayesInformed}
#' @return Returns \code{list} with prior and posterior samples (if applicable) from an object of class \code{bmult}
#' @export
samples.bmult <- function(x){
  output <- x$samples
  
  if(is.null(output)){
    cat("No prior or posterior samples were drawn.")
  }
  
  return(output)
}

#' S3 method for class 'bayes_factor.bmult'
#' 
#' @export
bayes_factor <- function (x) {
  UseMethod("bayes_factor")
}
#' Extracts information about computed Bayes factors from object of class \code{bmult}
#'
#' @param x object of class \code{bmult_bridge} as returned from \code{multBayesInformed}
#' @return Returns \code{list} with two dataframes from an object of class \code{bmult}. The first dataframe bf_all summarizes information
#' the Bayes factor for equality and inequality constraints. The second dataframe bf_ineq summarized information about the Bayes factor for inequality
#' constraints. 
#' @export
bayes_factor.bmult <- function(x){
  bf_list <- x$bf_list
  
  # Data frame 1
  bfer <- bf_list$logBFe_inequalities[,'logBFe_inequalities']
  bfe0 <- bf_list$logBFe_equalities[,'logBFe_equalities']
  bf_table <- data.frame(
    bf_type         = c('LogBFer', 'BFer', 'BFre'),
    bf_total        = as.numeric(bf_list$bf)
  )
  
  if(!purrr::is_empty(bfe0)){
    
    if(length(bfe0) == 1){
      
      bf_table$bf_equalities <- as.numeric(bf_list$equalities_list[[1]]$bf)
      
    } else {
      
      tab <- NULL
      
      for(i in 1:length(bfe0)){
        tab <- cbind(tab, as.numeric(bf_list$equalities_list[[i]]$bf))
      }
      
        colnames(tab) <- paste0('bf_eq_', 1:length(bfe0))
        bf_table      <- cbind(bf_table, tab)
        
    }
    
  }
  
  if(!purrr::is_empty(bfer)){
    
    if(length(bfer) == 1){
      
      bf_table$bf_inequalities <- c(bfer, exp(bfer), 1/exp(bfer))
      
    } else {
      
      tab <- NULL
      
      for(i in 1:length(bfer)){
        tab     <- cbind(tab, c(bfer[i], exp(bfer[i]), 1/exp(bfer[i])))
      }
      
      colnames(tab) <- paste0('bf_ineq_', 1:length(bfer))
      bf_table      <- cbind(bf_table, tab)
      
    }
    
  }
  
  row.names(bf_table) <- NULL
    
  # Data frame 2
  hypotheses <- sapply(x$restrictions$inequality_constraints$hyp, function(x) paste(x, collapse = ' '))
  bf_ineq_table <- data.frame(
    hyp = hypotheses,
    bf_list$logBFe_inequalities
  )
  
  output <- list(bf_table      = bf_table,
                 bf_ineq_table = bf_ineq_table)
  
  return(output)
  
}

#' print method for class \code{bmult_bridge}
#'
#' @param x object of class \code{bmult_bridge} as returned from \code{multBayesBfInequality}
#' @return The print methods print the results and return nothing
#' @export
print.bmult_bridge <- function(x){
  logml <- signif(x$logml,5)
  hyp   <- paste(x$hyp, collapse=" ")
  niter <- x$niter
  error <- x$error_measures$percentage

  output <- paste('Bridge sampling estimate of the log marginal likelihood:', logml,
               '\n\nHypothesis:\n\n', hyp,
               '\n\nEstimate obtained in', niter, 'iteration(s).',
               '\nPercentage Error:', error, sep = ' ')
  cat(output)
}

#' summary method for class \code{bmult_bridge}
#'
#' @param x object of class \code{bmult_bridge} as returned from \code{multBayesBfInequality} or \code{binomBayesBfInequality}
#' @return The summary method returns a \code{list} which contains the log marginal likelihood and associated error terms. 
#' @export
summary.bmult_bridge <- function(x){
  
  output <- list(logml   = x$logml,
                 hyp     = paste(x$hyp, collapse=" "),
                 re2     = x$error_measures$re2,
                 cv      = x$error_measures$cv,
                 percent = x$error_measures$percentage)
  class(output) <- c("summary.bmult_bridge", "list")
  
  printRes <- (paste('Bridge sampling log marginal likelihood estimate\n\n', signif(output$logml,5),
                     '\n\nHypothesis:\n\n', output$hyp,
                     '\n\nError Measures:\n\n',
                     'Relative Mean-Squared Error: ', round(output$re2,5),
                     '\nCoefficient of Variation: ', round(output$cv,5),
                     '\nPercentage Error: ', output$percent,
                     '\n\nNote:\nAll error measures are approximate.\n\n', sep=''))
  
  cat(printRes)
  invisible(output)
}

#' print method for class \code{bmult}
#'
#' @param x object of class \code{bmult} as returned from \code{multBayesInformed}
#' @return The print methods print the results and return nothing
#' @export
print.bmult <- function(x){
  
  hyp    <- paste(x$restrictions$full_model$hyp, collapse=" ")
  data   <- x$restrictions$full_model$counts_full
  
  factor_levels       <- x$restrictions$full_model$parameters_full
  a                   <- x$restrictions$full_model$alpha_full
  b                   <- x$restrictions$full_model$beta_full
  
  if(!is.null(b)){
    
    titleText <- 'Bayesian Evaluation of Order Constraints Between Independent Binomials'
    
  } else {
    
    titleText <- 'Bayesian Evaluation of Multinomial Order Constraints'
    
  }
  
  # Model specification
  hyp     <- paste(x$restrictions$full_model$hyp, collapse=" ")
  
  if(is.null(b)){
    
    prior <- data.frame(factor_levels=factor_levels, alpha = a)
    colnames(prior) <- c('', 'alpha')
    
  } else {
    
    prior <- data.frame(factor_levels=factor_levels, alpha = a, beta=b)
    colnames(prior) <- c('', 'alpha', 'beta')
    
  }
  
  
  output <- paste(titleText,
                  '\n\nHypothesis:\n\n', hyp,
                  '\n\nPrior Specification:\n\n')
  cat(output)
  print(prior)
}

#' summary method for class \code{bmult}
#'
#' @param x object of class \code{bmult} as returned from \code{multBayesInformed}
#' @return The summary method returns a \code{list} which contains the Bayes factor and associated hypotheses for the full
#' model, but also the separate for the independent equality and inequality constraints. 
#' @export
summary.bmult <- function(x){
  
  bf_list <- x$bf_list
  bf_type <- bf_list$bf_type
  bf      <- signif(bf_list$bf[bf_type], 5) 
  
  nr_equal      <- length(x$bf_list$logBFe_equalities[,'logBFe_equalities'])
  nr_inequal    <- length(x$bf_list$logBFe_inequalities[,'logBFe_inequalities'])
  eq_hyp_text   <- ifelse(nr_equal == 1, 'hypothesis', 'hypotheses')
  ineq_hyp_text <- ifelse(nr_inequal == 1, 'hypothesis.', 'hypotheses.')
  
  if(nr_equal == 0 & nr_inequal > 0){
    
    bfFootnote <- paste('\n\nBased on', nr_inequal, 'independent inequality-constrained', ineq_hyp_text)
    
  } else if(nr_equal > 0 &nr_inequal == 0){
    
    eq_hyp_text <- ifelse(nr_equal == 1, 'hypothesis.', 'hypotheses.')
    bfFootnote  <- paste('\n\nBased on', nr_equal, 'independent equality-constrained', eq_hyp_text)
    
  } else {
    
    bfFootnote <- paste('\n\nBased on', nr_equal, 'independent equality-constrained', eq_hyp_text, 'and',
                        nr_inequal, 'independent inequality-constrained', ineq_hyp_text)
    
  }
  
  bfText <- paste0('\n\nBayes factor estimate ', bf_type, ':\n\n', bf, bfFootnote)
  
  # 1. show model details
  hyp    <- paste(x$restrictions$full_model$hyp, collapse=" ")
  data   <- x$restrictions$full_model$counts_full
  total  <- x$restrictions$full_model$total_full
  
  # 2. show descriptives: if data is provided
  if(!is.null(data)){
    
    descriptivesText <- '\n\n2. Descriptives:\n'
    
    if(!is.null(total)){
      
      observed <- data.frame(factor_levels=factor_levels, 
                             observedCounts=data, 
                             totalN=total,
                             observedProportion=data/total)
      colnames(observed) <- c('', 'Counts', 'N', 'Proportions')
      
    } else {
      
      observed <- data.frame(factor_levels=factor_levels, observedCounts=data, observedProportion=data/sum(data))
      colnames(observed) <- c('', 'Counts', 'Proportions')
      
    }
    
  } 
  
  # 3. show prior or posterior estimates
  factor_levels       <- x$restrictions$full_model$parameters_full
  a                   <- x$restrictions$full_model$alpha_full
  b                   <- x$restrictions$full_model$beta_full
  
  if(!is.null(b)){
    
    estimates           <- .credibleIntervalPlusMedian(factor_levels=factor_levels, a=a, b=b, counts=counts, total=total)
    colnames(estimates) <- c('', 'alpha','beta', '2.5%', '50%', '97.5%')
    
  } else {
    
    estimates           <- .credibleIntervalPlusMedian(factor_levels=factor_levels, a=a, counts=counts)
    colnames(estimates) <- c('', 'alpha', '2.5%', '50%', '97.5%')
    
  }
  
  if(!is.null(data)){
    
    estimatesText <- '\n\n3. Posterior Median and Credible Intervals:\n'
    
  } else {
    
    estimatesText <- '\n\n3. Prior Median and Credible Intervals:\n'
    
  }
  
  output <- list(hyp = hyp, bf=bf, bf_type = bf_type, observed=observed, estimates=estimates)
  
  class(output) <- c("summary.bmult", "list")
  
  ## Print This ##
    printRes <- paste('Bayes factor analysis\n', 
                      '\nHypothesis:\n\n', hyp,
                      bfText, sep = ' ')
    
    cat(printRes)
    # print observed counts
    if(!is.null(data)){
      observed[,-1] <- signif(observed[,-1], 3)
      cat(descriptivesText)
      print(observed)
    }
    # print posterior estimates
    estimates[,-1] <- signif(estimates[,-1], 3)
    cat(estimatesText)
    print(estimates)
    
  ## Output ##
  invisible(output)
}

