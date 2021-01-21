#' @title S3 method for class \code{restriction_list.bmult}
#' 
#' @description Extracts restriction list from an object of class \code{bmult}
#' @param x  object of class \code{bmult} as returned from \code{\link{mult_bf_informed}} or \code{\link{binom_bf_informed}}
#' @param ... additional arguments passed to other methods
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
#' out_mult  <- mult_bf_informed(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
#' niter=1e3, seed=2020)
#' restriction_list <- restriction_list(out_mult)
#' @export
restriction_list <- function (x, ...) {
  UseMethod("restriction_list")
}

#' @title Extracts restriction list from an object of class \code{bmult}
#' @param restrictions specifies whether to extract restriction list for \code{equalities} or \code{inequalities}. Default is \code{inequalities}.
#' @param ... additional arguments, currently ignored
#' @inherit restriction_list
#' @export
restriction_list.bmult <- function(x, restrictions = 'inequalities', ...){
  
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
    
    message(paste0("\nNo restriction list created for ", expr, " constrained parameters.\n\n"))
    
  } 
  
  return(output)
}

#' @title S3 method for class \code{bridge_output.bmult}
#' 
#' @description Extracts bridge sampling output from object of class \code{bmult}
#' @param x object of class \code{bmult} as returned from \code{\link{mult_bf_informed}} or \code{\link{binom_bf_informed}}
#' @return Extracts output related to the bridge sampling routine. The output contains the following elements::
#' @return 
#' \describe{
#' \item{\code{$eval}}{
#' \itemize{
#' \item \code{q11}: log prior or posterior evaluations for prior or posterior samples
#' \item \code{q12}: log proposal evaluations for prior or posterior samples
#' \item \code{q21}: log prior or posterior evaluations for samples from proposal
#' \item \code{q22}: log proposal evaluations for samples from proposal
#' }}
#' \item{\code{$niter}}{number of iterations of the iterative updating scheme}
#' \item{\code{$logml}}{estimate of log marginal likelihood}
#' \item{\code{$hyp}}{evaluated inequality constrained hypothesis}
#' \item{\code{$error_measures}}{
#' \itemize{
#' \item \code{re2}: the approximate 
#' relative mean-squared error for the marginal likelihood estimate
#' \item \code{cv}: the approximate coefficient of variation for the marginal 
#' likelihood estimate (assumes that bridge estimate is unbiased)
#' \item \code{percentage}: the approximate percentage error of the marginal likelihood estimate
#' }}
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
#' out_mult  <- mult_bf_informed(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
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
    message("Bridge sampling was not applied.")
  } 
  
    return(output)
  
}

#' @title S3 method for class 'samples.bmult'
#' @description Extracts prior and posterior samples (if applicable) from an object of class \code{bmult}
#' @param x object of class \code{bmult} as returned from \code{\link{mult_bf_informed}} or \code{\link{binom_bf_informed}}
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
#' out_mult  <- mult_bf_informed(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
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
    message("No prior or posterior samples were drawn.")
  }
  
  return(output)
}

#' @title S3 method for class 'bayes_factor.bmult'
#' @description Extracts information about computed Bayes factors from object of class \code{bmult}
#' @param x object of class \code{bmult} as returned from \code{\link{mult_bf_informed}} or \code{\link{binom_bf_informed}}
#' or an object of class \code{bmult_bridge} as returned from \code{\link{mult_bf_inequality}} or \code{\link{binom_bf_inequality}}
#' @return Returns \code{list} with three \code{data.frames}. The first dataframe \code{bf_table} summarizes information
#' the Bayes factor for equality and inequality constraints. The second dataframe \code{error_measures} contains for the overall Bayes factor
#' the approximate relative mean-squared error \code{re2}, the approximate coefficient of variation \code{cv}, and the approximate percentage error \code{percentage}.
#' The third dataframe \code{$bf_ineq_table} summarized 
#' information about the Bayes factor for inequality constraints, that is, the log marginal likelihood estimates
#' for the constrained prior and posterior distribution. In addition, it contains for each independent Bayes factor
#' the approximate relative mean-squared error \code{re2}
#' 
#' @note 
#' In case the restricted hypothesis is tested against \eqn{H_0} four \code{data.frame}s will be returned.
#' The fourth dataframe \code{$bf_eq_table} summarizes information about the Bayes factor for equality constraints compared to the encompassing
#' hypothesis. 
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
#' out_mult  <- mult_bf_informed(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
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
  error_measures <- bf_list$error_measures
  
  # Data frame 1
  bfer <- bf_list$logBFe_inequalities[,'logBFe_inequalities']
  bfe0 <- bf_list$logBFe_equalities[,'logBFe_equalities']
  
  
  bf_table <- data.frame(
    bf_type         = colnames(bf_list$bf),
    bf_total        = as.numeric(bf_list$bf),
    stringsAsFactors = FALSE
  )
  
  # Separate output if alternative hypothesis is H_e or H_0
  # If alternative is H_0 then do not show breakdown
  if(colnames(bf_list$bf)[1] == "LogBFer"){
    
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
  }
  
  row.names(bf_table) <- NULL
  
  output <- list(bf_table        = bf_table,
                 error_measures  = error_measures)
    
  # Data frame 2
  if(!purrr::is_empty(bfer)){
    
    hypotheses   <- sapply(x$restrictions$inequality_constraints$hyp, function(x) paste(x, collapse = ' '))
    ineq_summary <- bf_list$logBFe_inequalities
    
    bf_ineq_table <- data.frame(
      hyp = hypotheses,
      ineq_summary,
      stringsAsFactors = FALSE
    )
    
    n_ineq <- nrow(ineq_summary)
    
    error_measures_prior <- error_measures_post <- NULL
    
    for(i in 1:n_ineq){
      
      bs_out <- x$bridge_output[[i]]
      error_measures_prior[i] <- ifelse(is.null(bs_out$prior), 0, bs_out$prior$error_measures$re2)
      error_measures_post[i]  <- ifelse(is.null(bs_out$post), 0, bs_out$post$error_measures$re2)
      
    }
    
    re2               <-  error_measures_prior + error_measures_post
    bf_ineq_table$re2 <- re2
    
    output$bf_ineq_table <- bf_ineq_table
    
  }
  
  # Include breakdown of equalities in separate list when alternative is H_0
  if(colnames(bf_list$bf)[1] == "LogBFr0"){
    
    if(!purrr::is_empty(bfe0)){
    
      hypotheses <- sapply(x$restrictions$equality_constraints$hyp, function(x) paste(x, collapse = ' '))
      eq_summary <- bf_list$logBFe_equalities
      
      bf_eq_table <- data.frame(
        hyp = hypotheses,
        eq_summary,
        stringsAsFactors = FALSE
      )
      
      output$bf_eq_table <- bf_eq_table
    
    }
    
  }
  
  return(output)
  
}

#' @title Print method for class \code{bmult_bridge}
#' @description Prints model specification
#'
#' @param x object of class \code{bmult_bridge} as returned from \code{\link{mult_bf_inequality}} 
#' or \code{\link{binom_bf_inequality}}
#' @param ... additional arguments, currently ignored
#' @return The print methods print the results from the bridge sampling
#' algorithm and return nothing
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
#' out_mult  <- mult_bf_inequality(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
#' niter=1e3, seed=2020)
#' out_mult
#' @export
print.bmult_bridge <- function(x, ...){
  logml <- signif(x$logml,5)
  hyp   <- paste(x$hyp, collapse=" ")
  niter <- x$niter
  error <- x$error_measures$percentage

  output <- paste('Bridge sampling estimate of the log marginal likelihood for\nthe constrained distribution:', logml,
               '\n\nHypothesis H_r:\n', hyp,
               '\n\nEstimate obtained in', niter, 'iteration(s).',
               '\nPercentage Error:', error,'\n', sep = ' ')
  cat(output)
}

#' @title summary method for class \code{bmult_bridge}
#' 
#' @description Summarizes bridge sampling results and associated error measures
#'
#' @param object object of class \code{bmult_bridge} as returned from \code{\link{mult_bf_inequality}} or \code{\link{binom_bf_inequality}}
#' @param ... additional arguments, currently ignored
#' @return  Invisibly returns a `list` which contains the log marginal likelihood and associated error terms. 
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
#' out_mult  <- mult_bf_inequality(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
#' niter=1e3, seed=2020)
#' summary(out_mult)
#' @export
summary.bmult_bridge <- function(object, ...){
  
  output <- list(logml   = object$logml,
                 hyp     = paste(object$hyp, collapse=" "),
                 re2     = object$error_measures$re2,
                 cv      = object$error_measures$cv,
                 percent = object$error_measures$percentage)
  class(output) <- c("summary.bmult_bridge", "list")
  
  return(output)
}

#' @title print method for class \code{summary.bmult_bridge}
#' 
#' @description Prints the summary of bridge sampling output
#'
#' @param x object of class \code{summary.bmult_bridge} as returned from \code{\link{summary.bmult_bridge}} 
#' @param ... additional arguments, currently ignored
#' @return The print methods print the summary of the bridge sampling output
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
#' out_mult  <- mult_bf_inequality(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
#' niter=1e3, seed=2020)
#' summary(out_mult)
#' @export
print.summary.bmult_bridge <- function(x, ...){
  
  printRes <- (paste('Bridge sampling log marginal likelihood estimate for \nthe constrained distribution: ', signif(x$logml,5),
                     '\n\nHypothesis H_r:\n', x$hyp,
                     '\n\nError Measures:\n\n',
                     'Relative Mean-Squared Error: ', round(x$re2,5),
                     '\nCoefficient of Variation: ', round(x$cv,5),
                     '\nPercentage Error: ', x$percent,
                     '\n\nNote:\nAll error measures are approximate.\n\n', sep=''))
  
  cat(printRes)
  
}

#' @title print method for class \code{bmult}
#' 
#' @description Prints model specification
#'
#' @param x object of class \code{bmult} as returned from \code{\link{mult_bf_informed}} or \code{\link{binom_bf_informed}}
#' @param ... additional arguments, currently ignored
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
#' out_binom  <- binom_bf_informed(x=x, n=n, Hr=Hr, a=a, b=b, niter=1e3,factor_levels, seed=2020)
#' out_binom
#' ## Multinomial Case
#' out_mult  <- mult_bf_informed(x=x, Hr=Hr, a=a, niter=1e3,factor_levels, seed=2020)
#' out_mult
#' @export
print.bmult <- function(x, ...){
  
  hyp    <- paste(x$restrictions$full_model$hyp, collapse=" ")
  counts <- x$restrictions$full_model$counts_full
  n      <- x$restrictions$full_model$total_full
  
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
  if(!is.null(counts)){
    
    descriptivesText <- '\n\n2. Descriptives:\n'
    
    if(!is.null(n)){
      
      observed <- data.frame(factor_levels=factor_levels, 
                             observedCounts=counts, 
                             totalN=n,
                             observedProportion=counts/n)
      colnames(observed) <- c('', 'Counts', 'N', 'Proportions')
      
    } else {
      
      observed <- data.frame(factor_levels=factor_levels, observedCounts=counts, 
                             observedProportion=counts/sum(counts))
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
  if(!is.null(counts)){
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
#' @param object object of class \code{bmult} as returned from \code{\link{mult_bf_informed}} or \code{\link{binom_bf_informed}}
#' @param ... additional arguments, currently ignored
#' @return   Invisibly returns a `list` which contains the Bayes factor and associated hypotheses for the full
#' model, but also the separate for the independent equality and inequality constraints. 
#' \describe{The summary method returns a \code{list} with the following elements:
#' \item{\code{$hyp}}{Vector containing the informed hypothesis as specified by the user}
#' \item{\code{$bf}}{Contains Bayes factor}
#' \item{\code{$bf}}{Contains relative mean-square error for the Bayes factor}
#' \item{\code{$bf_type}}{Contains Bayes factor type as specified by the user}
#' \item{\code{$cred_level}}{Credible interval for the posterior point estimates.}
#' \item{\code{$prior}}{List containing the prior parameters.}
#' \item{\code{$data}}{List containing the data.}
#' \item{\code{$nr_equal}}{Number of independent equality-constrained hypotheses.}
#' \item{\code{$nr_inequal}}{Number of independent inequality-constrained hypotheses.}
#' \item{\code{$estimates}}{Parameter estimates for the encompassing model
#' \itemize{
#' \item \code{factor_level}: Vector with category names
#' \item \code{alpha}: Vector with posterior concentration parameters of Dirichlet 
#' distribution (for multinomial models) or alpha parameters for independent beta 
#' distributions (for binomial models)
#' \item \code{beta}: Vector with beta parameters for independent beta 
#' distributions (for binomial models)
#' \item \code{lower}: Lower value of credible intervals of marginal beta distributions
#' \item \code{median}: Posterior median of marginal beta distributions
#' \item \code{upper}: Upper value of credible intervals of marginal beta distributions
#' }}
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
#' ## Binomial Case
#' out_binom  <- binom_bf_informed(x=x, n=n, Hr=Hr, a=a, b=b, niter=1e3,factor_levels, seed=2020)
#' summary(out_binom)
#' ## Multinomial Case
#' out_mult  <- mult_bf_informed(x=x, Hr=Hr, a=a, niter=1e3,factor_levels, seed=2020)
#' summary(out_mult)
#' @export
summary.bmult <- function(object, ...){
  
  bf_list <- object$bf_list
  bf_type <- bf_list$bf_type
  bf      <- signif(bf_list$bf[[bf_type]], 5) 
  re2     <- bf_list$error_measures$re2
  
  nr_equal      <- length(object$bf_list$logBFe_equalities[,'logBFe_equalities'])
  nr_inequal    <- length(object$bf_list$logBFe_inequalities[,'logBFe_inequalities'])
  
  # 1. show model details
  hyp    <- paste(object$restrictions$full_model$hyp, collapse=" ")
  counts <- object$restrictions$full_model$counts_full
  n      <- object$restrictions$full_model$total_full
  
  # 3. show prior or posterior estimates
  factor_levels       <- object$restrictions$full_model$parameters_full
  a                   <- object$restrictions$full_model$alpha_full
  b                   <- object$restrictions$full_model$beta_full
  cred_level          <- object$cred_level
  lower               <- ((1 - cred_level) / 2)
  upper               <- 1 - lower
  
  # for marginal beta distributions, determine b and total
  if(is.null(b))    b <- sum(a) - a
  if(!is.null(counts) & is.null(n)) n <- rep(sum(counts), length(counts))
  estimates_output <- .credibleIntervalPlusMedian(credibleIntervalInterval=cred_level, 
                                                  factor_levels=factor_levels, 
                                                  a=a, b=b, counts=counts, total=n)
  # slight changes for the output data frame
  colnames(estimates_output) <- c('factor_level', 'alpha','beta', 'lower', 'median', 'upper')
  estimates_output$alpha     <- ifelse(!is.null(counts), a + counts, a)
  estimates_output$beta      <- ifelse(!is.null(counts), b + n - counts, b) 
  
  output <- list(hyp = hyp, bf=bf, re2=re2,
                 bf_type = bf_type, 
                 cred_level=cred_level, 
                 prior=list(a=a, b=b), data=list(x=counts, n=n), 
                 nr_equal=nr_equal, 
                 nr_inequal=nr_inequal,
                 estimates=estimates_output)
  
  class(output) <- c("summary.bmult", "list")
    
  ## Output ##
  return(output)
}

#' @title print method for class \code{summary.bmult}
#' 
#' @description Prints the summary from Bayes factor analysis
#'
#' @param x object of class \code{bmult} or \code{summary.bmult} as returned from \code{\link{summary.bmult}} 
#' @param ... additional arguments, currently ignored
#' @return The print methods print the summary from the Bayes factor analysis and returns nothing
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
#' out_binom  <- binom_bf_informed(x=x, n=n, Hr=Hr, a=a, b=b, niter=1e3,factor_levels, seed=2020)
#' summary(out_binom)
#' ## Multinomial Case
#' out_mult  <- mult_bf_informed(x=x, Hr=Hr, a=a, niter=1e3,factor_levels, seed=2020)
#' summary(out_mult)
#' @export
print.summary.bmult <- function(x, ...){

  bf_type    <- x$bf_type
  bf         <- x$bf
  re2        <- x$re2
  nr_equal   <- x$nr_equal
  nr_inequal <- x$nr_inequal
  hyp        <- x$hyp

  a      <- x$prior$a
  b      <- x$prior$b
  counts <- x$data$x
  n      <- x$data$n
  factor_levels <- x$estimates$factor_level
  cred_level    <- x$cred_level

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
  re2Text <- paste0('\n\nRelative Mean-Square Error:\n\n', signif(re2, 3))

  # 3. show prior or posterior estimates
  lower               <- ((1 - cred_level) / 2)
  upper               <- 1 - lower
  lowerText           <- paste0(100*lower, '%')
  upperText           <- paste0(100*upper, '%')

  # for marginal beta distributions, determine b and total
  if(is.null(b))    b <- sum(a) - a
  if(!is.null(counts) & is.null(n)) n <- rep(sum(counts), length(counts))
  estimates <- .credibleIntervalPlusMedian(credibleIntervalInterval=cred_level,
                                           factor_levels=factor_levels,
                                           a=a, b=b, counts=counts, total=n)
  colnames(estimates) <- c('', 'alpha','beta', lowerText, '50%', upperText)

  if(!is.null(counts)){

    estimatesText <- '\n\nPosterior Median and Credible Intervals Of Marginal Beta Distributions:\n'

  } else {

    estimatesText <- '\n\nPrior Median and Credible Intervals Of Marginal Beta Distributions:\n'

  }

  ## Print This ##
  if(bf_type %in% c('BF0r', 'BFr0', 'LogBFr0')){

    printRes <- paste('Bayes factor analysis\n\n',
                      'Hypothesis H_0:\n',
                      'All parameters are exactly equal.\n\n',
                      'Hypothesis H_r:\n', hyp, bfText, re2Text, sep = ' ')

  } else {

    printRes <- paste('Bayes factor analysis\n\n',
                      'Hypothesis H_e:\n',
                      'All parameters are free to vary.\n\n',
                      'Hypothesis H_r:\n', hyp, bfText, re2Text, sep = ' ')

  }


  cat(printRes)

  # print posterior estimates
  estimates[,-c(1:3)] <- signif(estimates[,-c(1:3)], 3)
  cat(estimatesText)
  print(estimates)
}

#' Plot estimates
#' 
#' Plots the posterior estimates from the unconstrained multi- or binomial model.
#'
#' @param x A `summary.bmult`-object returned by `summary()`.
#' @param main `character`. A string used as title. Defaults to the informed
#'   hypothesis and the Bayes factor.
#' @param xlab `character`. A string used as title for the x-axis.
#' @param ylab `character`. A string used as title for the y-axis.
#' @param panel.first `expression`. Evaluated after the plot axes are set up but
#'   before any plotting takes place. This can be useful for drawing background
#'   grids or null distributions.
#' @inheritParams graphics::plot.window 
#' @param ... additional arguments passed to `plot.summary.bmult`-method.
#' @return Invisibly returns a `data.frame` with the plotted estimates.
#' @export
#'
#' @examples
#' # data
#' x <- c(3, 4, 10, 11, 7, 30)
#' # priors
#' a <- c(1, 1, 1, 1, 1, 1)
#' # restricted hypothesis
#' factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4', 'theta5', 
#'                    'theta6')
#' Hr <- c('theta1', '<',  'theta2', '&', 'theta3', '=', 
#' 'theta4', ',', 'theta5', '<', 'theta6')
#' output_total  <- mult_bf_informed(x, Hr, a, factor_levels, seed=2020, 
#' niter=1e3, bf_type = "BFer")
#' plot(summary(output_total))
#' 
#' # data for a big Bayes factor
#' x <- c(3, 4, 10, 11, 7, 30) * 1000
#' output_total  <- mult_bf_informed(x, Hr, a, factor_levels, seed=2020, 
#' niter=1e3, bf_type = "BFre")
#' plot(summary(output_total))

plot.summary.bmult <- function(x, main = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, panel.first = NULL, ...) {
  dat <- x$estimates
  # dat <- dat[order(dat$median, decreasing = FALSE), ]
  x_coord <- 1:length(dat$factor_level)
  
  # Set defaults
  if(is.null(xlab)) xlab <- "Parameter"
  if(is.null(ylab)) ylab <- bquote(widehat(theta)~"["~list(italic(q)[.(1-x$cred_level)], italic(q)[0.5], italic(q)[.(x$cred_level)])~"]")
  if(is.null(xlim)) xlim <- range(x_coord)
  if(is.null(ylim)) ylim <- range(dat[, c("upper", "lower")])
  
  op <- par(mar = c(4.25, 4.5, 2, 2) + 0.1)
  on.exit(expr = par(op))
  
  plot.new()
  plot.window(
    xlim = xlim
    , ylim = ylim
  )
  
  if(!is.null(panel.first)) panel.first
  
  axis(1, labels = dat$factor_level, at = x_coord)
  axis(2, las = 1)
  
  if(!is.null(main)) {
    title(main)
  } else {
    mtext(bquote(italic(H[r])~":"~.(x$hyp)))
    
    scientific_plotmath <- function(x) {
      x <- formatC(x, digits = 3, format = "g")
      parse(text = gsub("e\\+*(\\-*\\d+)", "%*%10^\\1", x))[[1]]
    }
    
    is_logbf <- grepl("Log", x$bf_type)
    
    mtext(
      bquote(.(if(is_logbf) "log(")*"BF"[italic(.(gsub("^.*BF", "", x$bf_type)))]*.(if(is_logbf) ")") == .(scientific_plotmath(x$bf)))
      , line = -1
    )
  }
  
  mtext(xlab, 1, cex = 1.4, line = 2.75)
  mtext(
    ylab
    , 2
    , cex = 1.3
    , line = 3
  )
  
  lines(
    x = x_coord
    , y = dat[, c("median")]
    , lwd = 1
    , col = grey(0.2)
  )
  
  points(
    x = x_coord
    , y = dat[, c("median")]
    , cex = 2.5
    , pch = 21
    , col = "white"
    , bg = "white"
  )
  
  arrows(
    x0 = x_coord
    , y0 = dat$lower
    , y1 = dat$upper
    , angle = 90
    , length = 0
    , code = 3
    , lwd = 1.25
  )
  
  points(
    x = x_coord
    , y = dat[, c("median")]
    , cex = 1.5
    , pch = 21
    , col = "black"
    , bg = "white"
    , lwd = 1.25
  )
  
  return(invisible(dat))
}

#' #' @export
#' #' @rdname plot.summary.bmult
#' 
#' plot.bmult <- function(x, ...) {
#'   plot(summary(x), ...)
#' }
