context("checks whether errors are identified correctly")

data(peas)
a             <- c(1, 1, 1, 1)     
counts        <- peas$counts
factor_levels <- peas$peas

test_that("incorrect specification of factor levels", {
  
  # incorrect factor name
  expect_error(multBayesInformed(factor_levels, c('roundYellow < notAFactor = roundGreen < wrinkledGreen'), a, counts),
               "\nThe following factor level(s) are invalid: notAFactor \n", fixed = TRUE)
  # using same factor name multiple times
  expect_error(multBayesInformed(factor_levels, c('roundYellow < wrinkledYellow = roundGreen < wrinkledGreen < roundYellow'), a, counts),
               "Do not use factor levels multiple times within the order restriction.", fixed = TRUE)
  # using opposing restriction signs
  expect_error(multBayesInformed(factor_levels, c('roundYellow < wrinkledYellow = roundGreen > wrinkledGreen'), a, counts),
               "Do not use the smaller and larger signs together within a restriction", fixed = TRUE)
})
