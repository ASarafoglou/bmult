context("checks whether categories get collapsed correctly")

data(peas)
a             <- c(1, 1, 1, 1)     
counts        <- peas$counts
factor_levels <- peas$peas

test_that("collapses categories correctly for Mendelian Peas example", {
  
  ## specify hypothesis explicitely
  Hr1                   <- c('roundYellow > wrinkledYellow = roundGreen > wrinkledGreen')
  output_total_explicit <- mult_bf_informed(factor_levels=factor_levels, 
                                             Hr=Hr1, a=a, x=counts,  
                                             niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(output_total_explicit$restrictions$equality_constraints$counts_equalities, c(101, 108))
  expect_equal(output_total_explicit$restrictions$inequality_constraints$counts_inequalities, 
               list(c(315, 209, 32)))
})

test_that("yields equal BF estimates for costraints with free parameters", {
  
  ## specify hypothesis explicitely
  Hr           <- c('roundYellow > wrinkledYellow , roundGreen > wrinkledGreen')
  output_total <- mult_bf_informed(factor_levels=factor_levels, 
                                    Hr=Hr, a=a, x=counts,  
                                    niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(output_total$restrictions$equality_constraints$counts_equalities, 
               numeric(0))
  expect_equal(output_total$restrictions$inequality_constraints$counts_inequalities, 
               list(c(315, 101, 108, 32)))
})


test_that("yields equal BF estimates for costraints with free and equal parameters", {
  
  ## specify hypothesis explicitely
  Hr           <- c('roundYellow > wrinkledYellow , roundGreen = wrinkledGreen')
  output_total <- mult_bf_informed(factor_levels=factor_levels, 
                                    Hr=Hr, a=a, x=counts,  
                                    niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(output_total$restrictions$equality_constraints$counts_equalities, c(108, 32))
  expect_equal(output_total$restrictions$inequality_constraints$counts_inequalities, 
               list(c(315, 101, 140)))
})

test_that("yields equal BF estimates for example with multiple equality constraints", {
  
  Hr            <- c('roundYellow = wrinkledYellow & roundGreen = wrinkledGreen')
  output_total  <- mult_bf_informed(factor_levels=factor_levels, 
                                     Hr=Hr, a=a, x=counts,  
                                     niter=5e3, bf_type = 'BFre', seed = 4)

  expect_equal(output_total$restrictions$inequality_constraints$counts_inequalities, NULL)
  expect_equal(output_total$restrictions$equality_constraints$counts_equalities, c(315, 101, 108, 32))

})

test_that("collapses categories correctly for ordered binomials", {
  
  ## specify hypothesis explicitly
  # priors
  a <- c(1, 1, 1, 1)
  b <- c(1, 1, 1, 1)
  x <- c(3, 7, 10, 12)
  n <- c(15, 12, 12, 12)
  # informed hypothesis
  factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4')
  Hr            <- c('theta1', '<',  'theta2', '=', 'theta3', '<', 'theta4')

  output_total_explicit <- binom_bf_informed(factor_levels=factor_levels, 
                                             Hr=Hr, a=a, b=b, x=x, n=n, 
                                             niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(output_total_explicit$restrictions$equality_constraints$counts_equalities, c(7, 10))
  expect_equal(output_total_explicit$restrictions$inequality_constraints$counts_inequalities, list(c(3, 17, 12)))
  
  expect_equal(
    output_total_explicit$bf_list,
    list(
      bf_type = "BFre",
      bf = structure(
        list(
          LogBFer = -1.78447390776466,
          BFer = 0.167885360956014,
          BFre = 5.95644548342722
        ),
        class = "data.frame",
        row.names = c(NA,-1L)
      ),
      error_measures = structure(
        list(
          re2 = 9.78994807280372e-06,
          cv = 0.00312888927141945,
          percentage = "0.3129%"
        ),
        class = "data.frame",
        row.names = c(NA,-1L)
      ),
      logBFe_equalities = structure(
        list(logBFe_equalities = -0.0207444369857086),
        row.names = c(NA,-1L),
        class = "data.frame"
      ),
      logBFe_inequalities = structure(
        list(
          logBFe_inequalities = -1.76372947077896,
          logml_prior = -1.79175946922805,
          logml_post = -0.0280299984490991
        ),
        class = "data.frame",
        row.names = c(NA,
                      -1L))))
})
