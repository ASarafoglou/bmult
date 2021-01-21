context("evaluate Bayes factors for equality constraints - Binomial")

test_that("yields equal BF estimates binomial equality constraints", {
  
  # Maarten Marsmans example
  factor_levels <- c('binom1', 'binom2', 'binom3', 'binom4')
  a             <- c(1, 1, 1, 1)
  b             <- c(1, 1, 1, 1)
  x             <- c(5, 10, 15, 14)
  n             <- c(17, 16, 16, 16)
  Hr            <- c('binom1', '=',  'binom2', '=', 'binom3', '=', 'binom4')
  output_total  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=x, n=n)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 6.16617565895481, 
                                                       BFer = 476.360851758647, 
                                                       BFre = 0.00209924891247499), class = "data.frame", row.names = c(NA, 
                                                                                                                                                 -1L)))
})

test_that("yields equal BFr0 estimates binomial equality constraints", {
  
  # Maarten Marsmans example
  factor_levels <- c('binom1', 'binom2', 'binom3', 'binom4')
  a             <- c(1, 1, 1, 1)
  b             <- c(1, 1, 1, 1)
  x             <- c(5, 10, 15, 14)
  n             <- c(17, 16, 16, 16)
  Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
  Hr1           <- c('binom1', '=',  'binom2', '=', 'binom3', '=', 'binom4')
  Hr2           <- c('binom1', '=',  'binom2')
  output_total  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, bf_type = 'BFr0',
                                     b=b, x=x, n=n, seed=2020)
  output_total2  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr1, a=a, bf_type = 'BFr0',
                                     b=b, x=x, n=n)
  output_total3  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr2, a=a, bf_type = 'BFr0',
                                     b=b, x=x, n=n)
  
  expect_equal(output_total$bf_list, list(bf_type = "BFr0", bf = structure(list(
    LogBFr0 = 8.04741180815038, BF0r = 0.000319928889188777, 
    BFr0 = 3125.69459586984), class = "data.frame", row.names = c(NA, 
                                                                  -1L)), bf0_table = structure(list(LogBFe0 = 6.16617565895481, 
                                                                                                    BFe0 = 476.360851758647, BF0e = 0.00209924891247499), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                              -1L)), bfr_table = structure(list(LogBFer = -1.88123614919558, 
                                                                                                                                                                                                                                BFer = 0.152401598156163, BFre = 6.56161098110872), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                        -1L)), error_measures = structure(list(re2 = 3.87905372316101e-05, 
                                                                                                                                                                                                                                                                                                                                                               cv = 0.00622820497668551, percentage = "0.6228%"), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                      -1L)), logBFe_inequalities = structure(list(logBFe_inequalities = -1.88123614919558, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  logml_prior = -3.17805383034795, logml_post = -1.29681768115237), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        -1L))))
  
  expect_equal(output_total2$bf_list, list(bf_type = "BFr0", bf = structure(list(
    LogBFr0 = 0, BF0r = 1, BFr0 = 1), class = "data.frame", row.names = c(NA, 
                                                                          -1L)), bf0_table = structure(list(LogBFe0 = 6.16617565895481, 
                                                                                                            BFe0 = 476.360851758647, BF0e = 0.00209924891247499), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                      -1L)), bfr_table = structure(list(LogBFer = 6.16617565895481, 
                                                                                                                                                                                                                                        BFer = 476.360851758647, BFre = 0.00209924891247499), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                  -1L)), error_measures = structure(list(re2 = 0, cv = 0, percentage = "0%"), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                  -1L)), logBFe_equalities = structure(list(logBFe_equalities = 6.16617565895481), row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 -1L), class = "data.frame")))
  
  expect_equal(output_total3$bf_list, list(bf_type = "BFr0", bf = structure(list(
    LogBFr0 = 0, BF0r = 1, BFr0 = 1), class = "data.frame", row.names = c(NA, 
                                                                          -1L)), bf0_table = structure(list(LogBFe0 = 0.843962315683342, 
                                                                                                            BFe0 = 2.32556336143926, BF0e = 0.430003334495739), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                    -1L)), bfr_table = structure(list(LogBFer = 0.843962315683342, 
                                                                                                                                                                                                                                      BFer = 2.32556336143926, BFre = 0.430003334495739), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                              -1L)), error_measures = structure(list(re2 = 0, cv = 0, percentage = "0%"), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                              -1L)), logBFe_equalities = structure(list(logBFe_equalities = 0.843962315683342), row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              -1L), class = "data.frame")))
  
  })

context("evaluate Bayes factors for equality constraints - Multinomial")

test_that("yields equal BF estimates multinomial equality constraints", {
  
  # Habermans Lifestresses
  data("lifestresses")
  a             <- rep(1, nrow(lifestresses))
  x             <- lifestresses$stress.freq
  output_total  <- mult_bf_equality(x, a)
  
  expect_equal(output_total$bf, structure(list(LogBFe0 = 3.29976435023366, 
                                               BFe0 = 27.1062505863656, BF0e = 0.0368918599351766), 
                                          class = "data.frame", row.names = c(NA, 
                                                                                                                                        -1L)))
})

test_that("yields equal BFr0 estimates multinomial equality constraints", {
  
  # Habermans Lifestresses
  data("peas")
  a             <- rep(1, nrow(peas))
  x             <- peas$counts
  Hr            <- c("1 < 2 < 3 < 4")
  Hr2           <- c("1 = 2 = 3 = 4")
  Hr3           <- c("1 = 2")
  
  output_total   <- mult_bf_informed(Hr=Hr, a=a, x=x, bf_type='BFr0', seed=2020)
  output_total2  <- mult_bf_informed(Hr=Hr2, a=a, x=x, bf_type='LogBFr0', seed=2020)
  output_total3  <- mult_bf_informed(Hr=Hr3, a=a, x=x, bf_type='LogBFr0', seed=2020)
  
  expect_equal(output_total$bf_list, list(bf_type = "BFr0", bf = structure(list(
    LogBFr0 = -15.3722077230551, BF0r = 4743129.57197044, BFr0 = 2.10831263372923e-07), class = "data.frame", row.names = c(NA, 
                                                                                                                            -1L)), bf0_table = structure(list(LogBFe0 = 142.870980852993, 
                                                                                                                                                              BFe0 = 1.11706542156224e+62, BF0e = 8.95202716597811e-63), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                             -1L)), bfr_table = structure(list(LogBFer = 158.243188576048, 
                                                                                                                                                                                                                                                                                               BFer = 5.29838603483751e+68, BFre = 1.88736719715189e-69), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                              -1L)), error_measures = structure(list(re2 = 9.86950050437589e-06, 
                                                                                                                                                                                                                                                                                                                                                                                                                                     cv = 0.00314157611787076, percentage = "0.3142%"), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            -1L)), logBFe_inequalities = structure(list(logBFe_inequalities = 158.243188576048, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        logml_prior = -3.17805383034795, logml_post = -161.421242406396), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              -1L))))
  
  expect_equal(output_total2$bf_list$bf, structure(list(LogBFr0 = 0, 
                                                        BF0r = 1, BFr0 = 1), class = "data.frame", row.names = c(NA, 
                                                                                                                 -1L)))
  
  expect_equal(output_total3$bf_list, list(bf_type = "LogBFr0", 
                                           bf = structure(list(LogBFr0 = 0, BF0r = 1, BFr0 = 1), class = "data.frame", row.names = c(NA, 
                                                                                                                                     -1L)), bf0_table = structure(list(LogBFe0 = 54.8269578439159, 
                                                                                                                                                                       BFe0 = 6.4721004805001e+23, BF0e = 1.54509344070432e-24), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                     -1L)), bfr_table = structure(list(LogBFer = 54.8269578439159, 
                                                                                                                                                                                                                                                                                                       BFer = 6.4721004805001e+23, BFre = 1.54509344070432e-24), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                     -1L)), error_measures = structure(list(re2 = 0, cv = 0, percentage = "0%"), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     -1L)), logBFe_equalities = structure(list(logBFe_equalities = 54.8269578439159), row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    -1L), class = "data.frame")))
  })
