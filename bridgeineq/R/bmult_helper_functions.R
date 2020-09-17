#' Split A Character Vector At A Given Symbol
#'
#' This function splits a character vector into separate components
#' using a specified breaking point.
#'
#' @param vector character vector to be broken into separate components.
#' @param signs character string containing a regular expression that is an element in vector and serves as breaking point.
#' @return a list with separate components
.splitAt        <- function(vector, signs){
  pos <- which(vector %in% signs)
  out <- unname(split(vector, cumsum(seq_along(vector) %in% pos)))
  return(out)
}

#' Collapse Categories That Are Constrained To Be Equal
#'
#' This function collapses categories that are constrained to be equal
#' It can collapse character vector (i.e., category names), as well as numeric vectors 
#' (i.e., data and concentration parameters). If the input is a character vector, the first element within the equalities 
#' is retained.
#'
#' @param counts either a character string or a numeric vector to be collapsed.
#' @param equalities index of elements in counts that are equality constrained.
#' @param is_numeric_value logical. If TRUE, the input is a numeric vector and if FALSE the input is a character vector.
#' @param correct logical. If TRUE a correction for marginalization is applied to numeric vectors. A correction is needed,
#' only if the function collapses concentration parameters.
#' @return collapsed character or numeric vector.
.collapseCategories <- function(counts, equalities, is_numeric_value = TRUE, correct = FALSE){
  
  if(is_numeric_value && !correct){
    new_counts <- lapply(equalities, function(x) c(rep(NA, length(x)-1), sum(counts[x])))
    for(i in seq_along(equalities)) counts[equalities[[i]]] <- new_counts[[i]]
    counts    <- counts[!is.na(counts)]
  } else if (is_numeric_value && correct) {
    new_counts <- lapply(equalities, function(x) c(rep(NA, length(x)-1), sum(counts[x]) - (length(x)-1)))
    for(i in seq_along(equalities)) counts[equalities[[i]]] <- new_counts[[i]]
    counts    <- counts[!is.na(counts)]
  } else {
    new_counts <- lapply(equalities, function(x) c(rep(NA, length(x)-1), counts[x[length(x)]]))
    for(i in seq_along(equalities)) counts[equalities[[i]]] <- new_counts[[i]]
    counts    <- counts[!is.na(counts)]
  }
  
  return(counts)
}

# convert a single string to OR
.splitString <- function(OR, factor_levels, signs = c(equal='=', smaller='<', larger='>', free=',', linebreak='&')){
  
  signs_c    <- stringr::str_c(signs, collapse = "")
  expression <- paste0('(?=[', signs_c, '])', sep='')
  
  # split order restriction at signs and trim away any whitespace
  OR_split <- strsplit(OR, expression, perl = TRUE)[[1]]
  
  return(OR_split)
}

# format hypotheses nicely 
.formatHypothesis <- function(hyp_list){
  
  hyp_list <- lapply(hyp_list, function(x) x[!x == '&'])
  hyp_vec  <- sapply(hyp_list, function(x) stringr::str_c(x, collapse = ' '))
  
  return(hyp_vec)
  
}