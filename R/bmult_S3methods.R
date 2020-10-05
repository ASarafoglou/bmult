#' @title S3 method for class \code{restriction_list.bmult}
#' 
#' @description Extracts restriction list from an object of class \code{bmult}
#' @param x object of class \code{bmult} as returned from \code{multBfInformed}
#' @param restrictions specifies whether to extract restriction list for \code{equalities} or \code{inequalities}. Default is \code{inequalities}.
#' @return Extracts restriction list and associated hypothesis from an object of class \code{bmult}
#' @examples 
#' # data
#' x <- c(3, 4, 10, 11)
#' n <- c(15, 12, 12, 12)
#' # priors
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' # informed hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4')
#' Hr            <- c('theta1', '<',  'theta2', '<', 'theta3', '<', 'theta4')
#' 
#' ## Multinomial Case
#' out_mult  <- multBfInformed(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
#' niter=1e3, seed=2020)
#' restriction_list <- restriction_list(out_mult)
#' @export
restriction_list <- function (x, restrictions = 'inequalities') {
  UseMethod("restriction_list")
}

#' @title Extracts restriction list from an object of class \code{bmult}
#'
#' @inherit restriction_list
#' @export
restriction_list.bmult <- function(x, restrictions = 'inequalities'){
  
  restrictions <- match.arg(restrictions, c('inequalities', 'equalities'))
                            
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

#' @title S3 method for class \code{bridge_output.bmult}
#' 
#' @description Extracts bridge sampling output from object of class \code{bmult}
#' @param x object of class \code{bmult_bridge} as returned from \code{multBfInformed}
#' @return Extracts restriction list from an object of class \code{bmult}. The bridge sampling output contains the following elements::
#' \describe{
#' \item{\code{$eval}}{list consisting of the following elements:
#' \itemize{
#' \item \code{q11}: log posterior evaluations for posterior samples
#' \item \code{q12}: log proposal evaluations for posterior samples
#' \item \code{q21}: log posterior evaluations for samples from proposal
#' \item \code{q22}: log proposal evaluations for samples from proposal
#' }}
#' \item{\code{$niter}}{number of iterations of the iterative updating scheme}
#' \item{\code{$logml}}{estimate of log marginal likelihood}
#' \item{\code{$hyp}}{character vector that contains the inequality constrained hypothesis }
#' \item{\code{$error_measures}}{error term of the bridge sampling estimate}
#' }
#' @examples 
#' # data
#' x <- c(3, 4, 10, 11)
#' n <- c(15, 12, 12, 12)
#' # priors
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' # informed hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4')
#' Hr            <- c('theta1', '<',  'theta2', '<', 'theta3', '<', 'theta4')
#' 
#' ## Multinomial Case
#' out_mult  <- multBfInformed(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
#' niter=1e3, seed=2020)
#' bridge_output <- bridge_output(out_mult)
#' @export
bridge_output <- function (x) {
UseMethod("bridge_output")
}

#' @title Extracts bridge sampling output from object of class \code{bmult}
#' @inherit restriction_list
#' @export
bridge_output.bmult <- function(x){
  output <- x$bridge_output
  
  if(is.null(output)){
    cat("Bridge sampling was not applied.")
  }
  
  return(output)
}

#' @title S3 method for class 'samples.bmult'
#' @description Extracts prior and posterior samples (if applicable) from an object of class \code{bmult}
#' @param x object of class \code{bmult_bridge} as returned from \code{multBfInformed}
#' @return Returns \code{list} with prior and posterior samples (if applicable) from an object of class \code{bmult}
#' @examples 
#' # data
#' x <- c(3, 4, 10, 11)
#' n <- c(15, 12, 12, 12)
#' # priors
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' # informed hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4')
#' Hr            <- c('theta1', '<',  'theta2', '<', 'theta3', '<', 'theta4')
#' 
#' ## Multinomial Case
#' out_mult  <- multBfInformed(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
#' niter=1e3, seed=2020)
#' sample_list <- samples(out_mult)
#' @export
samples <- function (x) {
  UseMethod("samples")
}

#' @title Extracts prior and posterior samples (if applicable) from an object of class \code{bmult}
#' 
#' @inherit samples
#' @export
samples.bmult <- function(x){
  output <- x$samples
  
  if(is.null(output)){
    cat("No prior or posterior samples were drawn.")
  }
  
  return(output)
}

#' @title S3 method for class 'bayes_factor.bmult'
#' @description Extracts information about computed Bayes factors from object of class \code{bmult}
#' @param x object of class \code{bmult_bridge} as returned from \code{multBfInformed}
#' @return Returns \code{list} with two \code{data.frames} from an object of class \code{bmult}. The first dataframe \code{bf_all} summarizes information
#' the Bayes factor for equality and inequality constraints. The second dataframe \code{bf_ineq} summarized information about the Bayes factor for inequality
#' constraints. 
#' @examples 
#' # data
#' x <- c(3, 4, 10, 11)
#' n <- c(15, 12, 12, 12)
#' # priors
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' # informed hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4')
#' Hr            <- c('theta1', '<',  'theta2', '<', 'theta3', '<', 'theta4')
#' 
#' ## Multinomial Case
#' out_mult  <- multBfInformed(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
#' niter=1e3, seed=2020)
#' bayes_factor(out_mult)
#' @export
bayes_factor <- function (x) {
  UseMethod("bayes_factor")
}
#' @title Extracts information about computed Bayes factors from object of class \code{bmult}
#'
#' @inherit bayes_factor
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
      
      bf_table$bf_equalities <- c(bfe0, exp(bfe0), 1/exp(bfe0))
      
    } else {
      
      tab <- NULL
      
      for(i in 1:length(bfe0)){
        tab <- cbind(tab, c(bfe0[i], exp(bfe0[i]), 1/exp(bfe0[i])))
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

#' @title Print method for class \code{bmult_bridge}
#' @description Prints model specification
#'
#' @param x object of class \code{bmult_bridge} as returned from \code{multBfBfInequality}
#' @return The print methods print the results and return nothing
#' @examples 
#' # data
#' x <- c(3, 4, 10, 11)
#' n <- c(15, 12, 12, 12)
#' # priors
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' # informed hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4')
#' Hr            <- c('theta1', '<',  'theta2', '<', 'theta3', '<', 'theta4')
#' 
#' ## Multinomial Case
#' out_mult  <- multBfInequality(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
#' niter=1e3, seed=2020)
#' out_mult
#' @export
print.bmult_bridge <- function(x){
  logml <- signif(x$logml,5)
  hyp   <- paste(x$hyp, collapse=" ")
  niter <- x$niter
  error <- x$error_measures$percentage

  output <- paste('Bridge sampling estimate of the log marginal likelihood for\nthe constrained distribution:', logml,
               '\n\nHypothesis H_r:\n', hyp,
               '\n\nEstimate obtained in', niter, 'iteration(s).',
               '\nPercentage Error:', error, sep = ' ')
  cat(output)
}

#' @title summary method for class \code{bmult_bridge}
#' 
#' @description Summarizes bridge sampling results and associated error measures
#'
#' @param x object of class \code{bmult_bridge} as returned from \code{multBfBfInequality} or \code{binomBayesBfInequality}
#' @return The summary method returns a \code{list} which contains the log marginal likelihood and associated error terms. 
#' @examples 
#' # data
#' x <- c(3, 4, 10, 11)
#' n <- c(15, 12, 12, 12)
#' # priors
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' # informed hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4')
#' Hr            <- c('theta1', '<',  'theta2', '<', 'theta3', '<', 'theta4')
#' 
#' ## Multinomial Case
#' out_mult  <- multBfInequality(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
#' niter=1e3, seed=2020)
#' summary(out_mult)
#' @export
summary.bmult_bridge <- function(x){
  
  output <- list(logml   = x$logml,
                 hyp     = paste(x$hyp, collapse=" "),
                 re2     = x$error_measures$re2,
                 cv      = x$error_measures$cv,
                 percent = x$error_measures$percentage)
  class(output) <- c("summary.bmult_bridge", "list")
  
  printRes <- (paste('Bridge sampling log marginal likelihood estimate for \nthe constrained distribution: ', signif(output$logml,5),
                     '\n\nHypothesis H_r:\n', output$hyp,
                     '\n\nError Measures:\n\n',
                     'Relative Mean-Squared Error: ', round(output$re2,5),
                     '\nCoefficient of Variation: ', round(output$cv,5),
                     '\nPercentage Error: ', output$percent,
                     '\n\nNote:\nAll error measures are approximate.\n\n', sep=''))
  
  cat(printRes)
  invisible(output)
}

#' @title print method for class \code{bmult}
#' 
#' @description Prints model specification
#'
#' @param x object of class \code{bmult} as returned from \code{multBfInformed}
#' @return The print methods print the model specifications and descriptives and return nothing
#' @examples 
#' # data
#' x <- c(3, 4, 10, 11)
#' n <- c(15, 12, 12, 12)
#' # priors
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' # informed hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4')
#' Hr            <- c('theta1', '<',  'theta2', '<', 'theta3', '<', 'theta4')
#' 
#' ## Binomial Case
#' out_binom  <- binomBfInformed(x=x, n=n, Hr=Hr, a=a, b=b, niter=1e3,factor_levels, seed=2020)
#' out_binom
#' ## Multinomial Case
#' out_mult  <- multBfInformed(x=x, Hr=Hr, a=a, niter=1e3,factor_levels, seed=2020)
#' out_mult
#' @export
print.bmult <- function(x){
  
  hyp    <- paste(x$restrictions$full_model$hyp, collapse=" ")
  data   <- x$restrictions$full_model$counts_full
  total  <- x$restrictions$full_model$total_full
  
  factor_levels       <- x$restrictions$full_model$parameters_full
  a                   <- x$restrictions$full_model$alpha_full
  b                   <- x$restrictions$full_model$beta_full
  
  # 1. show Model specification
  if(!is.null(b)){
    
    titleText <- 'Bayesian Evaluation of Order Constraints Between Independent Binomials'
    
  } else {
    
    titleText <- 'Bayesian Evaluation of Multinomial Order Constraints'
    
  }
  
  modelText <- paste(titleText, '\n\n1. Hypothesis:\n\n', hyp)
  
  # 2. show Descriptives: if data is provided
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
  
  # 3. show Prior specification
  hyp       <- paste(x$restrictions$full_model$hyp, collapse=" ")
  priorText <- '\n\n3. Prior Specification:\n\n'
  if(is.null(b)){
    
    prior <- data.frame(factor_levels=factor_levels, alpha = a)
    colnames(prior) <- c('', 'alpha')
    
  } else {
    
    prior <- data.frame(factor_levels=factor_levels, alpha = a, beta=b)
    colnames(prior) <- c('', 'alpha', 'beta')
    
  }
  
  # print model specification
  cat(modelText)
  # print observed counts
  if(!is.null(data)){
    observed[,-1] <- signif(observed[,-1], 3)
    cat(descriptivesText)
    print(observed)
  }
  # print prior
  cat(priorText)
  print(prior)
}

#' @title summary method for class \code{bmult}
#' 
#' @description Summarizes results from Bayes factor analysis
#'
#' @param x object of class \code{bmult} as returned from \code{multBfInformed}
#' @return The summary method returns a \code{list} which contains the Bayes factor and associated hypotheses for the full
#' model, but also the separate for the independent equality and inequality constraints. 
#' @examples 
#' # data
#' x <- c(3, 4, 10, 11)
#' n <- c(15, 12, 12, 12)
#' # priors
#' a <- c(1, 1, 1, 1)
#' b <- c(1, 1, 1, 1)
#' # informed hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4')
#' Hr            <- c('theta1', '<',  'theta2', '<', 'theta3', '<', 'theta4')
#' 
#' ## Binomial Case
#' out_binom  <- binomBfInformed(x=x, n=n, Hr=Hr, a=a, b=b, niter=1e3,factor_levels, seed=2020)
#' summary(out_binom)
#' ## Multinomial Case
#' out_mult  <- multBfInformed(x=x, Hr=Hr, a=a, niter=1e3,factor_levels, seed=2020)
#' summary(out_mult)
#' @export
summary.bmult <- function(x){
  
  bf_list <- x$bf_list
  bf_type <- bf_list$bf_type
  bf      <- signif(bf_list$bf[bf_type], 5) 
  
  nr_equal      <- length(x$bf_list$logBFe_equalities[,'logBFe_equalities'])
  nr_inequal    <- length(x$bf_list$logBFe_inequalities[,'logBFe_inequalities'])
  eq_hyp_text   <- ifelse(nr_equal == 1, 'hypothesis\n', 'hypotheses\n')
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
  
  # 3. show prior or posterior estimates
  factor_levels       <- x$restrictions$full_model$parameters_full
  a                   <- x$restrictions$full_model$alpha_full
  b                   <- x$restrictions$full_model$beta_full
  cred_level          <- x$cred_level
  lower               <- ((1 - cred_level) / 2)
  upper               <- 1 - lower
  lowerText           <- paste0(100*lower, '%')
  upperText           <- paste0(100*upper, '%')
  
  # for marginal beta distributions, determine b and total
  if(is.null(b))    b <- sum(a) - a
  if(!is.null(data) & is.null(total)) total <- rep(sum(data), length(data))
  estimates           <- .credibleIntervalPlusMedian(credibleIntervalInterval=cred_level, 
                                                     factor_levels=factor_levels, 
                                                     a=a, b=b, counts=data, total=total)
  colnames(estimates) <- c('', 'alpha','beta', lowerText, '50%', upperText)
  
  if(!is.null(data)){
    
    estimatesText <- '\n\nPosterior Median and Credible Intervals Of Marginal Beta Distributions:\n'
    
  } else {
    
    estimatesText <- '\n\nPrior Median and Credible Intervals Of Marginal Beta Distributions:\n'
    
  }
  
  output <- list(hyp = hyp, bf=bf, bf_type = bf_type, estimates=estimates)
  
  class(output) <- c("summary.bmult", "list")
  
  ## Print This ##
    printRes <- paste('Bayes factor analysis\n\n', 
                      'Hypothesis H_e:\n', 
                      'All parameters are free to vary.\n\n', 
                      'Hypothesis H_r:\n', hyp, bfText, sep = ' ')
    
    cat(printRes)

    # print posterior estimates
    estimates[,-c(1:3)] <- signif(estimates[,-c(1:3)], 3)
    cat(estimatesText)
    print(estimates)
    
  ## Output ##
  invisible(output)
}

