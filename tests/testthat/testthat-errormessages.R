context("checks whether errors are identified correctly")

data(peas)
a             <- c(1, 1, 1, 1)     
counts        <- peas$counts
factor_levels <- peas$peas

test_that("incorrect specification of factor levels", {
  
  # incorrect factor name
  expect_error(mult_bf_informed(factor_levels=factor_levels, Hr=c('roundYellow < notAFactor = roundGreen < wrinkledGreen'), a=a, x=counts),
               "\nThe following factor level(s) are invalid: notAFactor \n", fixed = TRUE)
  # using same factor name multiple times
  expect_error(mult_bf_informed(factor_levels=factor_levels, Hr=c('roundYellow < wrinkledYellow = roundGreen < wrinkledGreen < roundYellow'), a=a, x=counts),
               "Do not use factor levels multiple times within the order restriction.", fixed = TRUE)
  # using opposing restriction signs
  expect_error(mult_bf_informed(factor_levels=factor_levels, Hr=c('roundYellow < wrinkledYellow = roundGreen > wrinkledGreen'), a=a, x=counts),
               "Do not use the smaller and larger signs together within a restriction", fixed = TRUE)
})
