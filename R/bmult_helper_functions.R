# Split A Character Vector At A Given Symbol
.splitAt        <- function(vector, signs){
  pos <- which(vector %in% signs)
  out <- unname(split(vector, cumsum(seq_along(vector) %in% pos)))
  return(out)
}

# Collapse Categories That Are Constrained To Be Equal
# This function collapses categories that are constrained to be equal
# It can collapse character vector (i.e., category names), as well as numeric vectors 
# (i.e., data and concentration parameters). If the input is a character vector, the first element within the equalities 
# is retained.
.collapseCategories <- function(counts, equalities, is_numeric_value = TRUE, correct = FALSE, adjusted_priors_for_equalities=NULL){
  
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
.splitString <- function(OR, factor_levels, signs = c(equal='=', equal2='==', smaller='<', larger='>', free=',', linebreak='&')){
  
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

# compute credible interval plus median for beta distributions
.credibleIntervalPlusMedian <- function(credibleIntervalInterval = .95, factor_levels = NULL, a = NULL, b = NULL, counts = NULL, total = NULL) {
  
  if(credibleIntervalInterval > 0.99999 | credibleIntervalInterval < 0.00001 | 
     length(credibleIntervalInterval) != 1){
    stop("Credible interval must be a single number between 0 and 1.")
  }
  
  lower     <- (1 - credibleIntervalInterval) / 2
  upper     <- 1 - lower
  
  # for multinomial function
  if(is.null(b)){
    
    ciDf <- data.frame(factor_levels=factor_levels, 
                       alpha=as.character(a), 
                       lowerCI = NA, 
                       medianCI = NA, 
                       upperCI = NA)
    
    .checkAlphaAndData(alpha = a, counts = counts)
    
    if(is.null(counts)) {
      
      counts <- rep(0, length(a))
      
      } else {
      
      ciDf$alpha <- .str_align(a, counts)
      
      }
    
    total <- sum(counts)
    
    for(i in seq_along(a)){
      
      b <- sum(a) - a[i]
      
      binomResult <- stats::qbeta(c(lower, .5, upper), a[i] + counts[i] , b + total - counts[i])
      ciDf[i, -c(1:2)]   <- binomResult
      
    }
    
  } else {
    
    ciDf <- data.frame(factor_levels=factor_levels, 
                       alpha=as.character(a), 
                       beta=as.character(b), 
                       lowerCI = NA, 
                       medianCI = NA, 
                       upperCI = NA)
    
    .checkAlphaAndData(alpha = a, beta=b, counts = counts, total=total)
    
    for(i in seq_along(a)){
      
      if(is.null(counts)) {
        
        counts <- total <- rep(0, length(a))
        
      } else {
        
        ciDf$alpha <- .str_align(a, counts)
        ciDf$beta <- .str_align(b, total - counts)
          
      }
        
      
      binomResult <- stats::qbeta(c(lower, .5, upper), a[i] + counts[i] , b[i] + total[i] - counts[i])
      ciDf[i, -c(1:3)]   <- binomResult
      
    }
    
  }
  
  return(ciDf)
  
}

.str_align <- function(x, y) {
  aligned_x <- formatC(x, width = max(nchar(x)))
  aligned_y <- formatC(y, width = max(nchar(y)), flag = "-")
  paste(aligned_x, aligned_y, sep = ' + ')
}
