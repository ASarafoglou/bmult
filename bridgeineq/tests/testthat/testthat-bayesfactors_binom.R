context("evaluate Bayes factors for binomial variables")

factor_levels <- c('binom1', 'binom2', 'binom3', 'binom4')
a             <- c(1, 1, 1, 1)
b             <- c(1, 1, 1, 1)

test_that("yields equal BF estimates for Marsmans example", {
  
  # Maarten Marsmans example
  counts        <- c(3, 4, 10, 11)
  total         <- c(15, 12, 12, 12)
  Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
  output_total  <- binomBayesInformed(factor_levels, Hr, a, b, counts, total, niter = 5e4, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -2.56371376593299, 
                                                       BFer = 0.0770181811681734, BFre = 12.9839472295048), class = "data.frame", row.names = c(NA, 
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
  output_total  <- binomBayesInformed(factor_levels, Hr1, a=a, b=b, counts=counts, total=total, niter = 5e4, bf_type = 'LogBFer')
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 92.674128105046, 
                                                       BFer = 1.76954842211132e+40, BFre = 5.65115928733308e-41), class = "data.frame", row.names = c(NA, 
                                                                                                                                                      -1L)))
  
  output_total  <- binomBayesInformed(factor_levels, Hr2, a=a, b=b, counts=counts, total=total, niter = 5e4, bf_type = 'LogBFer')
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 52.0181386135618, 
                                                       BFer = 3.90113122080028e+22, BFre = 2.56335904485382e-23), class = "data.frame", row.names = c(NA, 
                                                                                                                                                      -1L)))
  
  output_total  <- binomBayesInformed(factor_levels, Hr3, a=a, b=b, counts=counts, total=total, niter = 5e4, bf_type = 'LogBFer')
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 366.160537081019, 
                                                       BFer = 1.05075325739535e+159, BFre = 9.51698215505738e-160), class = "data.frame", row.names = c(NA, 
                                                                                                                                                        -1L)))
  
})

test_that("BF estimate matches convergence expectations", {
  
  # yields equal BF estimates of 1 when no data is provided
  counts        <- c(0, 0, 0, 0)
  total         <- c(0, 0, 0, 0)
  Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
  output_total  <- binomBayesInformed(factor_levels, Hr, a, b, counts, total, niter = 5e4, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 0.00171697339499488, 
                                                       BFer = 1.00171844823778, BFre = 0.998284499760581), class = "data.frame", row.names = c(NA, 
                                                                                                                                               -1L)))
  
  # Very Large Values: BFre should converge to factorial(4)
  counts        <- c(3, 100, 800, 900)
  total         <- c(1e3, 1e3, 1e3, 1e3)
  Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
  output_total  <- binomBayesInformed(factor_levels, Hr, a, b, counts, total, niter = 5e4, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -3.16661950242359, 
                                                       BFer = 0.0421458312386302, BFre = 23.727139093259), class = "data.frame", row.names = c(NA, 
                                                                                                                                               -1L)))
  
  # Very Small Bayes factor
  Hr            <- c('binom1', '>',  'binom2', '>', 'binom3', '>', 'binom4')
  output_total  <- binomBayesInformed(factor_levels, Hr, a, b, counts, total, niter = 5e4, bf_type = 'LogBFer', seed=2020)
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 1589.69469924494, 
                                                       BFer = Inf, BFre = 0), class = "data.frame", row.names = c(NA, 
                                                                                                                  -1L)))
   
})

test_that("yields equal BF estimates for costraints with free parameters", {
  
  counts        <- c(3, 4, 10, 11)
  total         <- c(15, 12, 12, 12)
  Hr            <- c('binom1', '<',  'binom2', ',', 'binom3', ',', 'binom4')
  output_total  <- binomBayesInformed(factor_levels, Hr, a, b, counts, total, niter = 5e4, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -1.13707702062791, 
                                                       BFer = 0.320755213788234, BFre = 3.11764223000351), class = "data.frame", row.names = c(NA, 
                                                                                                                                               -1L)))
  
})