#' @title Prevalence of Statistical Reporting Errors
#'
#' @description This data set, "journals" provides a summary of statistical reporting errors (i.e., inconsistencies between reported test statistic
#' and reported p-value) of 16,695 research articles reporting results 
#' from null hypothesis significance testing (NHST). The selected articles were published in 
#' eight major journals in psychology between 1985 to 2013: 
#' \itemize{
#' \item *Developmental Psychology* (DP)
#' \item *Frontiers in Psychology* (FP)
#' \item *Journal of Applied Psychology* (JAP)
#' \item *Journal of Consulting and Clinical Psychology* (JCCP)
#' \item *Journal of Experimental Psychology: General* (JEPG)
#' \item *Journal of Personality and Social Psychology* (JPSP)
#' \item *Public Library of Science* (PLoS)
#' \item *Psychological Science* (PS)
#' }
#' 
#' In total, Nuijten et al. (2016) recomputed 258,105 p-values with the \code{R} software 
#' package \code{statcheck} which extracts statistics from articles and recomputes the p-values. 
#' The anonymized dataset and the data documentation was openly available on the 
#' Open Science Framework (\url{https://osf.io/d3ukb/}; \url{https://osf.io/c6ap2/}).
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
#' 
#'
#' @examples
#' data(journals)
#' # Prior specification 
#' # We assign a uniform Beta distribution on each binomial probability
#' a <- rep(1, 8)  
#' b <- rep(1, 8)  
#' 
#' x <- journals$errors 
#' n  <- journals$nr_NHST
#' factor_levels <- levels(journals$journal)
#' 
#' # restricted hypothesis
#' Hr1 <- c('JAP , PS , JCCP , PLOS , DP , FP , JEPG < JPSP')
#' out <- binom_bf_informed(x=x, n=n, Hr=Hr1, a=a, b=b, 
#' factor_levels=factor_levels, niter = 2e3)
#' 
#' summary(out)
'journals'
