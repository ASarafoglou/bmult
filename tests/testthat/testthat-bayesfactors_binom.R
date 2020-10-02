context("evaluate Bayes factors for informed hypotheses - Binomial")

factor_levels <- c('binom1', 'binom2', 'binom3', 'binom4')
a             <- c(1, 1, 1, 1)
b             <- c(1, 1, 1, 1)

test_that("yields equal BF estimates for Marsmans example", {
  
  # Maarten Marsmans example
  x <- c(5, 10, 15, 14)
  n <- c(17, 16, 16, 16)
  Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
  output_total  <- binomBfInformed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=x, n=n, niter = 5e3, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -1.88123614919558, 
                                                       BFer = 0.152401598156163, BFre = 6.56161098110872), class = "data.frame", row.names = c(NA, 
                                                                                                                                               -1L)))
})

test_that("yields equal BF estimates for Nuijten et al. 2016 example", {
  
  data(journals)
  a <- rep(1, 8)  
  b <- rep(1, 8)  
  counts <- journals$errors 
  total  <- journals$nr_NHST
  factor_levels <- levels(journals$journal)
  
  Hr1 <- c('JAP , PS , JCCP , PLOS , DP , FP , JEPG < JPSP')
  Hr2 <- c('JCCP < DP < JPSP')
  Hr3 <- c('JAP < PS < JCCP < PLOS < DP < FP < JEPG < JPSP')
  output_total  <- binomBfInformed(factor_levels=factor_levels, Hr=Hr1, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 133.707870599386, 
                                                       BFer = 1.17109031230873e+58, BFre = 8.53905108333245e-59), class = "data.frame", row.names = c(NA, 
                                                                                                                                                      -1L)))
  
  output_total  <- binomBfInformed(factor_levels=factor_levels, Hr=Hr2, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 63.5105935668673, 
                                                       BFer = 3.8220848883652e+27, BFre = 2.61637307701903e-28), class = "data.frame", row.names = c(NA, 
                                                                                                                                                     -1L)))
  
  output_total  <- binomBfInformed(factor_levels=factor_levels, Hr=Hr3, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 421.330436561272, 
                                                       BFer = 9.58260653962363e+182, BFre = 1.04355740357808e-183), class = "data.frame", row.names = c(NA, 
                                                                                                                                                        -1L)))
  
})

test_that("BF estimate matches convergence expectations", {
  
  # yields equal BF estimates of 1 when no data is provided
  counts        <- c(0, 0, 0, 0)
  total         <- c(0, 0, 0, 0)
  Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
  output_total  <- binomBfInformed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 0.00487829502509207, 
                                                       BFer = 1.00489021327864, BFre = 0.995133584530907), class = "data.frame", row.names = c(NA, 
                                                                                                                                               -1L)))
  # Very Large Values: BFre should converge to factorial(4)
  counts        <- c(3, 100, 800, 900)
  total         <- c(1e3, 1e3, 1e3, 1e3)
  Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
  output_total  <- binomBfInformed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -3.23721127565167, 
                                                       BFer = 0.0392732648369703, BFre = 25.4626144312463), class = "data.frame", row.names = c(NA, 
                                                                                                                                                -1L)))
  # Very Small Bayes factor
  Hr            <- c('binom1', '>',  'binom2', '>', 'binom3', '>', 'binom4')
  output_total  <- binomBfInformed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 1641.78003934897, 
                                                       BFer = Inf, BFre = 0), class = "data.frame", row.names = c(NA, 
                                                                                                                  -1L)))
})

test_that("yields equal BF estimates for costraints with free parameters", {
  
  counts        <- c(3, 4, 10, 11)
  total         <- c(15, 12, 12, 12)
  Hr            <- c('binom1', '<',  'binom2', ',', 'binom3', ',', 'binom4')
  output_total  <- binomBfInformed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -1.12851829614506, 
                                                       BFer = 0.323512250819858, BFre = 3.09107305045098), class = "data.frame", row.names = c(NA, 
                                                                                                                                               -1L)))
})