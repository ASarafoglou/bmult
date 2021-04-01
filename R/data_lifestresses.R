#' @title Memory of Life Stresses
#' 
#' @description This data set, "lifestresses", provides the number of reported life stresses (summed across participants) 
#' that occurred in specific months prior to an interview. This data set contains the subset of 147 participants who 
#' reported one negative life event over the time span of 18 months prior to an interview. Description taken from the JASP (2020)
#' data library.
#'
#' @docType data
#'
#' @usage data(lifestresses)
#'
#' @format A \code{data.frame} with 18 rows and 3 variables:
#' \describe{
#'   \item{\code{month}}{The month in which participants reported a stressful life event.}
#'   \item{\code{stress.freq}}{The number of participants who reported a life stress in the particular month prior to an interview.}
#'   \item{\code{stress.percentage}}{The percentage of participants who reported a life stress in the particular month prior to an interview.}
#' }
#'
#' @keywords datasets
#'
#' @references 
#' \insertRef{haberman1978}{multibridge}
#' 
#' \insertRef{jasp}{multibridge}
#' 
#' \insertRef{sarafoglou2020evaluatingPreprint}{multibridge}
#' 
#' \insertRef{uhlenhuth1977remembering}{multibridge}
#'
#' @examples
#' data(lifestresses)
#' # Prior specification 
#' # We assign a uniform Dirichlet distribution, that is, we set all 
#' # concentration parameters to 1
#' a             <- rep(1, 18)
#' x             <- lifestresses$stress.freq
#' factor_levels <- lifestresses$month
#' # Test the following restricted Hypothesis:
#' # Hr: month1 > month2 > ... > month18 
#' Hr            <- paste0(1:18, collapse=">"); Hr
#' out  <- mult_bf_informed(x=x, Hr=Hr, a=a, factor_levels=factor_levels,
#' niter=1e3, bf_type = 'BFre', seed = 4)
#' m1 <- summary(out)
#' m1
'lifestresses'