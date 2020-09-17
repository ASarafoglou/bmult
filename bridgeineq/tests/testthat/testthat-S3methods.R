context("evaluate S3 methods")

test_that("yields equal Bayes factor output for Haberman example", {
  
  data(lifestresses)
  a             <- rep(1, 18)
  counts        <- lifestresses$stress.freq
  factor_levels <- lifestresses$month
  Hr            <- paste0(1:18, collapse=">")
  output_total  <- multBayesInformed(factor_levels, Hr, a, counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(bayes_factor(output_total),
               list(
                 bf_table = structure(
                   list(
                     bf_type = structure(
                       c(3L, 1L, 2L),
                       .Label = c("BFer", "BFre",
                                  "LogBFer"),
                       class = "factor"
                     ),
                     bf_total = c(-5.14219588563597,
                                  0.00584484100893979, 171.091052514599),
                     bf_inequalities = c(-5.14219588563597,
                                         0.00584484100893979, 171.091052514599)
                   ),
                   row.names = c(NA,-3L),
                   class = "data.frame"
                 ),
                 bf_ineq_table = structure(
                   list(
                     hyp = structure(1L, .Label = "1 > 2 > 3 > 4 > 5 > 6 > 7 > 8 > 9 > 10 > 11 > 12 > 13 > 14 > 15 > 16 > 17 > 18", class = "factor"),
                     logBFe_inequalities = -5.14219588563597,
                     logml_prior = -36.3954452080331,
                     logml_post = -31.2532493223971
                   ),
                   class = "data.frame",
                   row.names = c(NA,-1L)
                 )
               ))
  
})

test_that("yields equal Bayes factor output for Mendelian Peas example", {
  
  data(peas)
  a             <- c(1, 1, 1, 1)     
  counts        <- peas$counts
  factor_levels <- levels(peas$peas)
  Hr            <- c('roundYellow > wrinkledYellow = roundGreen > wrinkledGreen')
  output_total  <- multBayesInformed(factor_levels, Hr, a, counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(bayes_factor(output_total),
               list(
                 bf_table = structure(
                   list(
                     bf_type = structure(
                       c(3L, 1L, 2L),
                       .Label = c("BFer", "BFre",
                                  "LogBFer"),
                       class = "factor"
                     ),
                     bf_total = c(-4.07162056514532,
                                  0.0170497358993748, 58.6519348981041),
                     bf_equalities = c(-2.33226548597986,
                                       0.0970755744318043, 10.3012524608083),
                     bf_inequalities = c(-1.73935507916546,
                                         0.175633633889566, 5.69367027177025)
                   ),
                   row.names = c(NA,-3L),
                   class = "data.frame"
                 ),
                 bf_ineq_table = structure(
                   list(
                     hyp = structure(1L, .Label = "roundYellow > wrinkledYellow = roundGreen > wrinkledGreen", class = "factor"),
                     logBFe_inequalities = -1.73935507916546,
                     logml_prior = -1.79187528451288,
                     logml_post = -0.0525202053474192
                   ),
                   class = "data.frame",
                   row.names = c(NA,-1L)
                 )
               ))
  
})

test_that("yields equal Bayes factor output for costraints with free parameters", {
  
  data(peas)
  a             <- c(1, 1, 1, 1)     
  counts        <- peas$counts
  factor_levels <- levels(peas$peas)
  
  Hr           <- c('roundYellow > wrinkledYellow , roundGreen > wrinkledGreen')
  output_total <- multBayesInformed(factor_levels, Hr, a, counts,  niter=5e3, bf_type = 'BFre', seed = 4)

  expect_equal(bayes_factor(output_total),
               list(
                 bf_table = structure(
                   list(
                     bf_type = structure(
                       c(3L, 1L, 2L),
                       .Label = c("BFer", "BFre",
                                  "LogBFer"),
                       class = "factor"
                     ),
                     bf_total = c(-2.40735611810453,
                                  0.0900530697897938, 11.1045631463119),
                     bf_inequalities = c(-2.40735611810453,
                                         0.0900530697897938, 11.1045631463119)
                   ),
                   row.names = c(NA,-3L),
                   class = "data.frame"
                 ),
                 bf_ineq_table = structure(
                   list(
                     hyp = structure(1L, .Label = "roundYellow > wrinkledYellow , roundGreen > wrinkledGreen", class = "factor"),
                     logBFe_inequalities = -2.40735611810453,
                     logml_prior = -2.48180647354589,
                     logml_post = -0.0744503554413648
                   ),
                   class = "data.frame",
                   row.names = c(NA,-1L)
                 )
               ))
  
})

test_that("tests S3 methods for example with multiple equality and inequality constraints", {
  
  data(lifestresses)
  a             <- rep(1, 18)
  counts        <- lifestresses$stress.freq
  factor_levels <- lifestresses$month
  Hr            <- c('1 < 2 & 5 < 3 < 4 < 8 & 10 = 11 = 12 & 13 = 14 = 15 < 16')
  output_total  <- multBayesInformed(factor_levels, Hr, a, counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(bayes_factor(output_total),
               list(
                 bf_table = structure(
                   list(
                     bf_type = structure(
                       c(3L, 1L, 2L),
                       .Label = c("BFer", "BFre", "LogBFer"),
                       class = "factor"
                     ),
                     bf_total = c(2.58119285238096,13.2128897901512, 0.0756836707095968),
                     bf_eq_1 = c(-1.03504256068524, 0.355211262568491, 2.81522605102416),
                     bf_eq_2 = c(-2.20127503702833, 0.110661970261721, 9.03652806501591),
                     bf_ineq_1 = c(0.816349968384871, 2.2622275475873, 0.442042181418272),
                     bf_ineq_2 = c(5.03993011644694, 154.459220489584, 0.00647420074263186),
                     bf_ineq_3 = c(-0.0387696347372691, 0.961972288622288, 1.03953098423674)
                   ),
                   class = "data.frame",
                   row.names = c(NA,-3L)
                 ),
                 bf_ineq_table = structure(
                   list(
                     hyp = structure(
                       3:1,
                       .Label = c("& 13 = 14 = 15 < 16",
                                  "& 5 < 3 < 4 < 8", "1 < 2"),
                       class = "factor"
                     ),
                     logBFe_inequalities = c(0.816349968384871, 5.03993011644694,-0.0387696347372691),
                     logml_prior = c(-0.693147180559945,-3.17805383034795,-0.287587820248094),
                     logml_post = c(-1.50949714894482,-8.21798394679488,-0.248818185510824)
                   ),
                   class = "data.frame",
                   row.names = c(NA,-3L)
                 )
               ))
})