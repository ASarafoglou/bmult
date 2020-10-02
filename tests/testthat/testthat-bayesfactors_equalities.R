context("evaluate Bayes factors for equality constraints - Binomial")

test_that("yields equal BF estimates binomial equality constraints", {
  
  # Maarten Marsmans example
  factor_levels <- c('binom1', 'binom2', 'binom3', 'binom4')
  a             <- c(1, 1, 1, 1)
  b             <- c(1, 1, 1, 1)
  x             <- c(5, 10, 15, 14)
  n             <- c(17, 16, 16, 16)
  Hr            <- c('binom1', '=',  'binom2', '=', 'binom3', '=', 'binom4')
  output_total  <- binomBfInformed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=x, n=n)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 6.16617565895481, 
                                                       BFer = 476.360851758647, 
                                                       BFre = 0.00209924891247499), class = "data.frame", row.names = c(NA, 
                                                                                                                                                 -1L)))
})

context("evaluate Bayes factors for equality constraints - Multinomial")

test_that("yields equal BF estimates multinomial equality constraints", {
  
  # Habermans Lifestresses
  data("lifestresses")
  a             <- rep(1, nrow(lifestresses))
  x             <- lifestresses$stress.freq
  output_total  <- multBfEquality(x, a)
  
  expect_equal(output_total$bf, structure(list(LogBFe0 = 3.29976435023366, 
                                               BFe0 = 27.1062505863656, BF0e = 0.0368918599351766), 
                                          class = "data.frame", row.names = c(NA, 
                                                                                                                                        -1L)))
})