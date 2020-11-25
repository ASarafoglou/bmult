context("evaluate S3 method bayes_factor")

test_that("yields equal Bayes factor output for Haberman example", {
  
  data(lifestresses)
  a             <- rep(1, 18)
  counts        <- lifestresses$stress.freq
  factor_levels <- lifestresses$month
  Hr            <- paste0(1:18, collapse=">")
  output_total  <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                     x=counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  
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
  output_total  <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                     x=counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  
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
  output_total <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                    x=counts,  niter=5e3, bf_type = 'BFre', seed = 4)

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

test_that("tests Bayes factor output for example with multiple equality and inequality constraints", {
  
  data(lifestresses)
  a             <- rep(1, 18)
  counts        <- lifestresses$stress.freq
  factor_levels <- lifestresses$month
  Hr            <- c('1 < 2 & 5 < 3 < 4 < 8 & 10 = 11 = 12 & 13 = 14 = 15 < 16')
  output_total  <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                     x=counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(bayes_factor(output_total), list(bf_table = structure(list(
    bf_type = structure(c(3L, 1L, 2L), .Label = c("BFer", "BFre", "LogBFer"), class = "factor"), 
    bf_total = c(7.39408166027864,  1626.33071352535, 0.000614881088872957), 
    bf_eq_1 = c(-2.20127503702833, 0.110661970261721, 9.03652806501591), 
    bf_eq_2 = c(0.0675301291989676, 1.06986249321224, 0.934699558442806), 
    bf_ineq_1 = c(0.816349968384871, 2.2622275475873, 0.442042181418272), 
    bf_ineq_2 = c(4.78908069462129, 120.190825787053, 0.00832010258230312), 
    bf_ineq_3 = c(3.92239590510185, 50.5213442348234, 0.0197936142663187)), 
    class = "data.frame", row.names = c(NA, -3L)), 
    bf_ineq_table = structure(list(hyp = structure(3:1, .Label = c("& 13 = 14 = 15 < 16", 
                                                                   "& 5 < 3 < 4 < 8", "1 < 2"), class = "factor"), 
                                   logBFe_inequalities = c(0.816349968384871, 4.78908069462129, 3.92239590510185), 
                                   logml_prior = c(-0.693147180559945, -3.17805383034795, -0.287587820248094), 
                                   logml_post = c(-1.50949714894482, -7.96713452496923, -4.20998372534994)), 
                              class = "data.frame", row.names = c(NA, -3L)))) 
  })
