context("evaluate restriction list")

test_that("yields equal restriction list for Haberman example", {
  
  data(lifestresses)
  a             <- rep(1, 18)
  counts        <- lifestresses$stress.freq
  factor_levels <- lifestresses$month
  Hr            <- paste0(1:18, collapse=">")
  output_total  <- mult_bf_informed(factor_levels=factor_levels, 
                                     Hr=Hr, a=a, x=counts,  
                                     niter=5e2, bf_type = 'BFre', seed = 4)
  
  expect_equal(output_total$restrictions$inequality_constraints$boundaries,
               list(
                 list(
                   list(lower = 2:18, upper = integer(0)),
                   list(lower = 3:18,
                        upper = 1L),
                   list(lower = 4:18, upper = 1:2),
                   list(lower = 5:18,
                        upper = 1:3),
                   list(lower = 6:18, upper = 1:4),
                   list(lower = 7:18,
                        upper = 1:5),
                   list(lower = 8:18, upper = 1:6),
                   list(lower = 9:18,
                        upper = 1:7),
                   list(lower = 10:18, upper = 1:8),
                   list(lower = 11:18, upper = 1:9),
                   list(lower = 12:18, upper = 1:10),
                   list(lower = 13:18, upper = 1:11),
                   list(lower = 14:18,
                        upper = 1:12),
                   list(lower = 15:18, upper = 1:13),
                   list(lower = 16:18, upper = 1:14),
                   list(lower = 17:18,
                        upper = 1:15),
                   list(lower = 18L, upper = 1:16),
                   list(lower = integer(0), upper = 1:17)
                 )
               ))
})

test_that("yields equal restriction list for Mendelian Peas example", {
  
  data(peas)
  a             <- c(1, 1, 1, 1)     
  counts        <- peas$counts
  factor_levels <- levels(peas$peas)
  Hr            <- c('roundYellow > wrinkledYellow = roundGreen > wrinkledGreen')
  output_total  <- mult_bf_informed(factor_levels=factor_levels, 
                                     Hr=Hr, a=a, x=counts,  
                                     niter=5e2, bf_type = 'BFre', seed = 4)
  
  expect_equal(output_total$restrictions$inequality_constraints$boundaries,
               list(list(
                 list(lower = 2:3, upper = integer(0)),
                 list(lower = 3L,
                      upper = 1L),
                 list(lower = integer(0), upper = 1:2)
               )))
})

test_that("yields equal restriction list when factor levels are mixed up", {
  
  x <- c(1, 4, 1, 10)
  a <- c(1, 1, 1, 1)
  factor_levels <- c('mult1', 'mult2', 'mult3', 'mult4')
  Hr <- c('mult2 > mult1 , mult3 = mult4')
  restrictions <- generate_restriction_list(x=x, Hr=Hr, a=a, factor_levels=factor_levels)
  
  expect_equal(restrictions$full_model, 
               structure(list(hyp = c("mult2", ">", "mult1", ",", "mult3", "=", "mult4"), 
                              parameters_full = c("mult2", "mult1", "mult3", "mult4"), 
                              alpha_full = c(1, 1, 1, 1), counts_full = c(4, 1, 1, 10)), 
                         class = c("list", "bmult_rl_full")))
  
  expect_equal(restrictions$inequality_constraints$counts_inequalities, 
               list(c(4, 1, 11)))
})