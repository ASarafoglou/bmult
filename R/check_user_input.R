# Checks user input and makes sure that analyses can be executed correctly
.checkRestrictionListClass <- function(restrictions){
  
  # If restrictions are not a restriction list, check whether order restriction is specified correctly
  if(!inherits(restrictions, 'bmult_rl') & !inherits(restrictions, 'bmult_rl_ineq')){
    
    # stop function if OR is not a character vector
    if(!is.character(restrictions)) stop('This analysis requires the specification of restrictions either as character vector or as object of class bmult_rl as returned from generateRestrictionList.')
    
  }
    
} 
.checkAlphaAndData <- function(alpha, beta=NULL, counts=NULL, total=NULL){
  
  # stop function if alpha is empty or values are not numeric
  if(is.null(alpha)) stop('Please specify alpha parameters.')
  if(!is.numeric(alpha)) stop('Alpha needs to be a numeric vector.')
  if(any(alpha < 0)) stop('Alpha cannot be negative.')
  
  if(!is.null(counts)){
    
    # stop function if data values are not numeric
    if(!is.numeric(counts)) stop('Data needs to be a numeric vector.')
    if(any(counts < 0)) stop('Counts cannot be negative.')
    # stop function if alpha and data are not of the same length
    if(length(alpha) != length(counts)) stop('alpha and counts are not of the same length.')
    
  }
  
  if(!is.null(total)){
    # stop function if total and counts are not of the same length
    if(length(total) != length(counts)) stop('counts and total are not of the same length.')
    if(any(total < counts)) stop('Total number of observations cannot be smaller than number of successes.')
  }
  
  if(!is.null(beta)){
    if(is.null(counts) & !is.null(total) || !is.null(counts) & is.null(total)) stop('For ordered binomials, please specify number of successes and total number of observations.')
  }
  
}
.checksIfMatrix <- function(mat){
  
  # stops function is provided input is not a matrix
  if(!is.matrix(mat)) stop('The bridge sampling method requires a matrix of samples.')
  
}
.checkNrParameters <- function(samples = NULL, boundaries = NULL, factors_analysis = NULL, alpha = NULL, counts = NULL){
  
  compare <- NULL
  
  if(!is.null(samples)) compare <- c(compare, ncol(samples))
  if(!is.null(boundaries)) compare <- c(compare, length(boundaries))
  if(!is.null(factors_analysis)) compare <- c(compare, length(factors_analysis))
  if(!is.null(alpha)) compare <- c(compare, length(alpha))
  if(!is.null(counts)) compare <- c(compare, length(counts))
  
  if(length(unique(compare)) != 1) {
    
    stop('The number of parameters differ.')
    
  }
  
}
.checkOrderRestriction <- function(OR, signs = c(equal='=', smaller='<', larger='>', free=',', linebreak='&')){
  # check whether OR is specified
  anyRestrictionPresent <- any(signs %in% OR)
  if(!anyRestrictionPresent){
    # check whether user input was numeric representation
    stop(paste0('No valied order restriction found. Valid restriction signs are:', stringr::str_c(signs, collapse= ' '), '\n'))
  }
  
  # checks if smaller and larger sign were used within the same restriction  
  distinct_restrictions   <- .splitAt(OR, signs['linebreak']) 
  check_smaller_larger    <- distinct_restrictions %>% purrr::map(function(x) signs['smaller'] %in% x & signs['larger'] %in% x) 
  smaller_larger_violated <- any(unlist(check_smaller_larger))
  
  if(smaller_larger_violated){
    
    stop('Do not use the smaller and larger signs together within a restriction')

  }
  
}
.checkFactorLevelsInOR <- function(OR, factor_levels, signs = c(equal='=', smaller='<', larger='>', free=',', linebreak='&')){
  # check whether factor levels are present in OR
  anyFactorLevelsPresent <- any(factor_levels %in% OR)
  if(!anyFactorLevelsPresent){
    # check whether user input was numeric representation
    stop(paste('No match between order restrictions and factor levels.'))
  }
  
  # checks for invalid factor levels in OR
  OR_nosigns         <- purrr::discard(OR, function(x) any(signs %in% x))                 # discard signs from the order restriction
  OR_invalid_factors <- purrr::discard(OR_nosigns, function(x) any(factor_levels %in% x)) # identifies invalid factor levels
  
  if(!purrr::is_empty(OR_invalid_factors)){
    stop(paste('\nThe following factor level(s) are invalid:', stringr::str_c(OR_invalid_factors, collapse= ' '), '\n'))
  }
  
  # stop if factor levels were used multiple times
  if(length(OR_nosigns) != length(unique(OR_nosigns))){
    
    stop('Do not use factor levels multiple times within the order restriction.')
    
  }
  
}
.checkForNumericRepresentation <- function(OR, factor_levels, signs = c(equal='=', smaller='<', larger='>', free=',', linebreak='&')){
  # returns order restriction of error message
  
  OR_no_signs          <- purrr::discard(OR, function(x) any(signs %in% x))                  # discard signs
  OR_invalid_factors   <- purrr::discard(OR_no_signs, function(x) any(factor_levels %in% x)) # identifies invalid factor levels
  
  if(!purrr::is_empty(OR_invalid_factors)){
    
    numbers_only <- all(!grepl("\\D", OR_invalid_factors))
    
    if(numbers_only){
      
      OR_indeces              <- as.numeric(OR_invalid_factors)
      OR_corresponds_to_index <- OR_indeces %in% seq_along(factor_levels)
      
      if(!all(OR_corresponds_to_index)){
        
        invalid_indeces <- stringr::str_c(OR_invalid_factors[!OR_corresponds_to_index], collapse = ' ')
        stop(paste('\nThe following indeces are invalid:', invalid_indeces, '\n'))
        
      } else {
        
        # if numbers only, convert to regular OR
        position_in_OR     <- match(OR_invalid_factors, OR)
        OR[position_in_OR] <- factor_levels[OR_indeces]
        
      }
    } 
  }
  return(OR)
}
.checkSpecifiedConstraints <- function(OR, factor_levels, signs = c(equal='=', smaller='<', larger='>', free=',', linebreak='&')){
  # returns order restriction or error message
  
  # if necessary, split string into character vector
  if(length(OR)==1) OR <- .splitString(OR, factor_levels, signs)
  
  # if necessary, transform index representation to factor level representation
  # trim OR of whitespace, before proceeding
  OR <- stringr::str_trim(OR)
  if(any(OR == '')) OR <- OR[-which(OR == '')]
  OR <- .checkForNumericRepresentation(OR, factor_levels, signs)
  
  # check factor levels in order restriction (do they exist?)
  .checkFactorLevelsInOR(OR, factor_levels, signs)
  
  # check the order restriction
  .checkOrderRestriction(OR, signs)
  
  return(OR)
}

# .checkRestrictionSigns <- function(signs, sings_default = c(equal='=', smaller='<', larger='>', free=c',', linebreak='&')){
#   # returns restriction signs or error message
#   
#   # stop if same symbol was used multiple times
#   if(length(signs) != length(unique(signs))){
#     
#     stop('Do not match multiple character strings to the same restrictions.')
#     
#   }
#   
#   # stop if restriction symbols are numbers
#   numbers_as_restrictions <- any(!grepl("\\D", signs))
#   if(numbers_as_restrictions){
#     
#     stop('Numbers cannot be used as restriction signs.')
#     
#   }
#   
#   # stop if input is not a named vector
#   names_correct <- c('equal', 'smaller', 'larger', 'free', 'linebreak')
#   nr_signs      <- length(signs)
#   sign_names    <- names(signs)
#   nr_names      <- length(sign_names)
#   
#   if(is.null(sign_names) | !(all(sign_names %in% names_correct))){
#     
#     stop('Please provide a named vector that matches character strings to the following restrictions: equal, smaller, larger, free, linebreak.')
#     
#   }
#   
#   # correct if not all components in symbol are specified
#   if(nr_signs < 5) {
#     
#     signs_default[sign_names] <- signs
#     signs                     <- signs_default 
#     
#   }
#   
#   return(signs)
#   
# }
