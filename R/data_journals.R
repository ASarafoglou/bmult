#' Prevalence of Statistical Reporting Errors
#'
#' This data set, "journals" provides a summary of statistical reporting errors (i.e., inconsistencies between reported test statistic
#' and reported p-value) of 16,695 research articles reporting results 
#' from null hypothesis significance testing (NHST). The selected articles were published in 
#' eight major journals in psychology between 1985 to 2013: *Developmental Psychology* (DP), the *Frontiers in Psychology* (FP), the 
#' *Journal of Applied Psychology* (JAP), the *Journal of Consulting and Clinical Psychology* (JCCP), 
#' *Journal of Experimental Psychology: General* (JEPG), the *Journal of Personality and Social Psychology* (JPSP), 
#' the *Public Library of Science* (PLoS), *Psychological Science* (PS). 
#' 
#' In total, \insertCite{nuijten2016prevalence; textual}{multibridge} recomputed 258,105 p-values with the software package `statcheck` 
#' which extracts statistics from articles and recomputes the p-values \insertCite{epskamp2016statcheck}{multibridge}. 
#' The anonymized dataset and the data documentation was openly available on the Open Science Framwework (\url{https://osf.io/d3ukb/}; 
#' \url{https://osf.io/c6ap2/}).
#'
#' @docType data
#'
#' @usage data(journals)
#'
#' @format A \code{data.frame} with 8 rows and 14 variables:
#' 
#' |Variable Name | Description|
#' | --- | --- |
#' | `journal` | The journal name a research article was published in.|
#' | `articles_downloaded` | The number of articles downloaded per journal.|
#' | `articles_with_NHST` | The number of articles with NHST results.|
#' | `perc_articles_with_NHST` | The percentage of all downloaded articles that had NHST results.|
#' | `nr_NHST` | The total number of NHST results.|
#' | `mean_nr_NHST_per_article_with_NHST` | The mean number of NHST results per article that had at least one NHST result.|
#' | `mean_nr_NHST_per_article_all_included` | The mean number of NHST results in all downloaded articles.|
#' | `errors` | The total number of errors.|
#' | `dec_errors` | The total number of decision errors (i.e., an error that may have changed the statistical conclusion of the result).|
#' | `perc_errors` | The percentage of all results that was an error.|
#' | `perc_dec_errors` | The percentage of all results that was a decision error.|
#' | `perc_articles_with_errors` | The percentage of all articles that had at least one error.|
#' | `perc_articles_with_dec_errors` | The percentage of all articles that had at least one error.|
#' | `APAfactor` | APA factor: number of detected NHST results / total number of detected p values.|
#' 
#'
#' @keywords datasets
#'
#' @references 
#' \insertRef{nuijten2016prevalence}{multibridge} 
#' \insertRef{epskamp2016statcheck}{multibridge}
#'
#'
#' @examples
#' data(journals)
#' # Prior specification 
#' # We assign a uniform Beta distribution on each binomial probability
#' a <- rep(1, 8)  
#' b <- rep(1, 8)  
#' 
#' counts <- journals$errors 
#' total  <- journals$nr_NHST
#' factor_levels <- levels(journals$journal)
#' 
#' # restricted hypothesis
#' Hr1 <- c('JAP , PS , JCCP , PLOS , DP , FP , JEPG < JPSP')
#' Hr2 <- c('JCCP < DP < JPSP')
#' Hr3 <- c('JAP < PS < JCCP < PLOS < DP < FP < JEPG < JPSP')
#' 
#' out <- binomBayesInformed(factor_levels, Hr1, a=a, b=b, counts=counts, total=total, niter = 5e4, bf_type = 'LogBFer'); out
#' out <- binomBayesInformed(factor_levels, Hr2, a=a, b=b, counts=counts, total=total, niter = 5e4, bf_type = 'LogBFer'); out
#' out <- binomBayesInformed(factor_levels, Hr3, a=a, b=b, counts=counts, total=total, niter = 5e4, bf_type = 'LogBFer'); out
#' 
#' summary(out)
'journals'
