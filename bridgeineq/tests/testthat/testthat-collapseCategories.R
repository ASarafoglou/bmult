context("checks whether categories get collapsed correctly")

data(peas)
a             <- c(1, 1, 1, 1)     
counts        <- peas$counts
factor_levels <- peas$peas

test_that("collapses categories correctly for Mendelian Peas example", {
  
  ## specify hypothesis explicitely
  Hr1                   <- c('roundYellow > wrinkledYellow = roundGreen > wrinkledGreen')
  output_total_explicit <- multBayesInformed(factor_levels, Hr1, a, counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(output_total_explicit$restrictions$equality_constraints$counts_equalities, c(101, 108))
  expect_equal(output_total_explicit$restrictions$inequality_constraints$counts_inequalities, 
               list(c(315, 209, 32)))
})

test_that("yields equal BF estimates for costraints with free parameters", {
  
  ## specify hypothesis explicitely
  Hr           <- c('roundYellow > wrinkledYellow , roundGreen > wrinkledGreen')
  output_total <- multBayesInformed(factor_levels, Hr, a, counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(output_total$restrictions$equality_constraints$counts_equalities, 
               numeric(0))
  expect_equal(output_total$restrictions$inequality_constraints$counts_inequalities, 
               list(c(315, 101, 108, 32)))
})


test_that("yields equal BF estimates for costraints with free and equal parameters", {
  
  ## specify hypothesis explicitely
  Hr           <- c('roundYellow > wrinkledYellow , roundGreen = wrinkledGreen')
  output_total <- multBayesInformed(factor_levels, Hr, a, counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(output_total$restrictions$equality_constraints$counts_equalities, c(108, 32))
  expect_equal(output_total$restrictions$inequality_constraints$counts_inequalities, 
               list(c(315, 101, 140)))
})

test_that("yields equal BF estimates for example with multiple equality constraints", {
  
  Hr            <- c('roundYellow = wrinkledYellow & roundGreen = wrinkledGreen')
  output_total  <- multBayesInformed(factor_levels, Hr, a, counts,  niter=5e3, bf_type = 'BFre', seed = 4)

  expect_equal(output_total$restrictions$inequality_constraints$counts_inequalities, NULL)
  expect_equal(output_total$restrictions$equality_constraints$counts_equalities, c(315, 101, 108, 32))

})