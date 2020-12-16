#' @title Mendelian Laws of Inheritance
#'
#' @description This data set, "peas", provides the categorization of crossbreeds between a plant variety that produced 
#' round yellow peas with a plant variety that produced wrinkled green peas. This data set contains the categorization
#' of 556 plants that were categorized either as (1) round and yellow, (2) wrinkled and yellow, (3) round and green, 
#' or (4) wrinkled and green.
#'
#' @docType data
#'
#' @usage data(peas)
#'
#' @format A \code{data.frame} with 4 rows and 2 variables:
#' \describe{
#'   \item{\code{peas}}{Crossbreeds that are categorized as 'roundYellow', 'wrinkledYellow', 'roundGreen', or 'wrinkledGreen'.}
#'   \item{\code{counts}}{The number of plants assigned to a one of the crossbreed categories.}
#' }
#'
#' @keywords datasets
#'
#' @references 
#' \insertRef{mulder2020generalizationPreprint}{multibridge}
#' 
#' \insertRef{robertson1978testing}{multibridge}
#' 
#' \insertRef{sarafoglou2020evaluatingPreprint}{multibridge}
#'
#'
#' @examples
#' data("peas")
#' # Prior specification 
#' # We assign a uniform Dirichlet distribution, that is, we set all 
#' # concentration parameters to 1
#' a <- c(1, 1, 1, 1)     
#' 
#' x <- peas$counts
#' factor_levels <- levels(peas$peas)
#' # Test the following mixed Hypothesis:
#' # Hr: roundYellow > wrinkledYellow = roundGreen > wrinkledGreen 
#' #
#' # Be careful: Factor levels are usually ordered alphabetically!
#' # When specifying hypotheses using indexes, make sure they refer to the 
#' # correct factor levels.
#' Hr <- c('1 > 2 = 3 > 4') 
#' # To avoid mistakes, write out factor levels explicitly:
#' Hr <- c('roundYellow > wrinkledYellow = roundGreen > wrinkledGreen')
#' 
#' out <- mult_bf_informed(x=x, Hr=Hr, a=a, factor_levels=factor_levels, niter=2e3,
#' bf_type = 'BFre')
#' summary(out)
'peas'