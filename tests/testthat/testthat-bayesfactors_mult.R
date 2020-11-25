context("evaluate Bayes factors for informed hypotheses - Multinomial")

test_that("yields equal BF estimates of 1 when no data is provided", {
  
  a             <- c(1, 1, 1, 1, 1)
  counts        <- c(0, 0, 0, 0, 0)
  factor_levels <- paste0('g', 1:5)
  Hr            <- c('g1 < g2 = g3 = g4 < g5')
  output_total  <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                     x=counts,  niter=5e3, bf_type = 'BFre', 
                                     seed = 4)
  
  expect_equal(output_total$bf_list$bf, 
               structure(list(LogBFer = 0, BFer = 1, BFre = 1),
                         class = "data.frame", row.names = c(NA, -1L)))
  
})

test_that("yields equal BF estimates for Haberman example", {
  
  data(lifestresses)
  a             <- rep(1, 18)
  counts        <- lifestresses$stress.freq
  factor_levels <- lifestresses$month
  Hr            <- paste0(1:18, collapse=">")
  output_total  <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                     x=counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  
  expect_equal(output_total$bf_list$bf, 
               data.frame(LogBFer = -5.14220, BFer=0.00584, BFre=171.09105),
               tolerance = 0.002)
  
})

test_that("yields equal BF estimates for Mendelian Peas example", {
  
  data(peas)
  a             <- c(1, 1, 1, 1)     
  counts        <- peas$counts
  factor_levels <- peas$peas
  
  ## specify hypothesis explicitely
    Hr1                   <- c('roundYellow > wrinkledYellow = roundGreen > wrinkledGreen')
    output_total_explicit <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr1, a=a, 
                                               x=counts,  niter=5e3, bf_type = 'BFre', seed = 4)
    # total BF
    expect_equal(output_total_explicit$bf_list$bf,
                 data.frame(LogBFer=-4.071621, BFer=0.01704974, BFre=58.65193),
                 tolerance = 0.002)
    # inequality BF
    expect_equal(output_total_explicit$bf_list$logBFe_inequalities,
                 data.frame(logBFe_inequalities=-1.73936, logml_prior=-1.79188, logml_post=-0.05252),
                 tolerance = 0.002)
    # equality BF
    expect_equal(output_total_explicit$bf_list$logBFe_equalities,
                 data.frame(logBFe_equalities=-2.33227),
                 tolerance = 0.002)
  
  ## specify hypothesis using indeces
    Hr2                <- c('1 > 2 = 3 > 4')
    output_total_index <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr2, a=a, 
                                            x=counts,  niter=5e3, bf_type = 'BFre', seed = 4)
    # total BF
    expect_equal(output_total_index$bf_list$bf,
                 data.frame(LogBFer=-4.071621, BFer=0.01704974, BFre=58.65193),
                 tolerance = 0.002)
    # inequality BF
    expect_equal(output_total_index$bf_list$logBFe_inequalities,
                 data.frame(logBFe_inequalities=-1.73936, logml_prior=-1.79188, logml_post=-0.05252),
                 tolerance = 0.002)
    # equality BF
    expect_equal(output_total_index$bf_list$logBFe_equalities,
                 data.frame(logBFe_equalities=-2.33227),
                 tolerance = 0.002)
    
    ## specify hypothesis using mix of explicit names and indeces
    Hr3              <- c('roundYellow > 2 = 3 > wrinkledGreen')
    output_total_mix <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr3, a=a, 
                                          x=counts,  niter=5e3, bf_type = 'BFre', seed = 4)
    # total BF
    expect_equal(output_total_mix$bf_list$bf,
                 data.frame(LogBFer=-4.071621, BFer=0.01704974, BFre=58.65193),
                 tolerance = 0.002)
    # inequality BF
    expect_equal(output_total_mix$bf_list$logBFe_inequalities,
                 data.frame(logBFe_inequalities=-1.73936, logml_prior=-1.79188, logml_post=-0.05252),
                 tolerance = 0.002)
    # equality BF
    expect_equal(output_total_mix$bf_list$logBFe_equalities,
                 data.frame(logBFe_equalities=-2.33227),
                 tolerance = 0.002)
})

test_that("yields equal BF estimates for costraints with free parameters", {
  
  data(peas)
  a             <- c(1, 1, 1, 1)     
  counts        <- peas$counts
  factor_levels <- levels(peas$peas)
  
  Hr           <- c('roundYellow > wrinkledYellow , roundGreen > wrinkledGreen')
  output_total <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                    x=counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  # total BF
  expect_equal(output_total$bf_list$bf,
               data.frame(LogBFer=-2.407356, BFer=0.09005307, BFre=11.10456),
               tolerance = 0.002)
  # inequality BF
  expect_equal(output_total$bf_list$logBFe_inequalities,
               data.frame(logBFe_inequalities=-2.407356, logml_prior=-2.481806, logml_post=-0.07445036),
               tolerance = 0.002)
})

test_that("yields equal BF estimates for costraints with free and equal parameters", {
  
  data(peas)
  a             <- c(1, 1, 1, 1)     
  counts        <- peas$counts
  factor_levels <- levels(peas$peas)
  
  Hr           <- c('roundYellow > wrinkledYellow , roundGreen = wrinkledGreen')
  output_total <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                    x=counts, niter=5e3, bf_type = 'BFre', seed = 4)
  # total BF
  expect_equal(output_total$bf_list$bf,
               data.frame(LogBFer=18.19921, BFer=80133886, BFre=1.247912e-08),
               tolerance = 0.002)
  # inequality BF
  expect_equal(output_total$bf_list$logBFe_inequalities,
               data.frame(logBFe_inequalities=-1.16125, logml_prior=-0.8742136, logml_post= 0.2870361),
               tolerance = 0.002)
  # equality BF
  expect_equal(output_total$bf_list$logBFe_equalities,
               data.frame(logBFe_equalities=19.36046),
               tolerance = 0.002)
})

test_that("yields equal BF estimates for example with multiple equality constraints", {
  
  data(peas)
  a             <- c(1, 1, 1, 1)     
  counts        <- peas$counts
  factor_levels <- levels(peas$peas)
  Hr            <- c('roundYellow = wrinkledYellow & roundGreen = wrinkledGreen')
  output_total  <- mult_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                     x=counts,  niter=5e3, bf_type = 'BFre', seed = 4)
  # equality BF
  expect_equal(output_total$bf_list$logBFe_equalities, structure(list(
    logBFe_equalities = c(54.8269578439159, 19.36045897615)), 
    row.names = c(NA, -2L), class = "data.frame"))
})

