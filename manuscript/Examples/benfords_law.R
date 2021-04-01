library('multibridge')

# Observed counts
x <- c(509, 353, 177, 114,  77,  77,  53,  73,  64)
# Concentration parameters
a <-  rep(1, 9)

# Execute the analysis

# Equality constrained hypotheses
 
  # Expected proportions for H_0 and H_r2
  p0  <- log10((1:9 + 1)/1:9)
  pr2 <- rep(1/9, 9)
  
  results_H0_He   <- mult_bf_equality(x = x, a = a, p = p0)
  results_Hr2_He  <- mult_bf_equality(x = x, a = a, p = pr2)
  
  logBFe0  <- results_H0_He$bf$LogBFe0
  logBFer2 <- results_Hr2_He$bf$LogBFe0
  
# Inequality constrained hypotheses
  
  # Labels for categories of interest
  factor_levels <- 1:9
  # Specifying the informed Hypothesis
  Hr1 <- c('1 > 2 > 3 > 4 > 5 > 6 > 7 > 8 > 9')
  Hr3 <- c('1 > 2 = 3 = 4 = 5 = 6 = 7 > 8 > 9')
  
  results_He_Hr1 <- mult_bf_informed(x = x, Hr = Hr1, a = a, 
                                     factor_levels = factor_levels, 
                                     bf_type = 'LogBFer', seed = 2020)
  results_He_Hr3 <- mult_bf_informed(x = x, Hr = Hr3, a = a, 
                                     factor_levels = factor_levels, 
                                     bf_type = 'LogBFer', seed = 2020)
  
  logBFer1 <- summary(results_He_Hr1)$bf
  logBFer3 <- summary(results_He_Hr3)$bf
  
# Compute Bayes Factors and Posterior Model Probabilites
  bayes_factors <- data.frame(
    BFType = c('LogBFe0', 'LogBFr10', 'LogBFr20', 'LogBFr30'), 
    LogBF  = c(logBFe0, -logBFer1 + logBFe0, -logBFer2 + logBFe0, -logBFer3 + logBFe0)
    )
  
  denominator <- sum(exp(-logBFe0), 1, exp(-logBFer1), exp(-logBFer2), exp(-logBFer3))
  post_probs <- data.frame(
    Hyps = c('p(H0 | x)', 'p(He | x)', 'p(Hr1 | x)', 'p(Hr2 | x)', 'p(Hr3 | x)'), 
    Prob = c(exp(-logBFe0), 1, exp(-logBFer1), exp(-logBFer2), exp(-logBFer3))/denominator
  )
  
# Create Plot
  first_digits <- 1:9
  benford <- log10((first_digits + 1) / first_digits)
  
  plot(
    summary(results_He_Hr1) # could be substituted with: summary(results_He_Hr3)
    , xlab = "Leading digit"
    , main = ""
    , panel.first = {
      lines(x = first_digits, y = benford, lty = "22", col = grey(0.5))
      points(x = first_digits, y = benford, pch = 16, col = "white", cex = 2)
      points(x = first_digits, y = benford, pch = 16, bg = "white", col = grey(0.7), cex = 0.8)
    }
  )
  
  points(x = first_digits, y = benford, pch = 16, bg = "white", col = grey(0.7), cex = 0.8)
  
  legend(
    "right"
    , legend = c("Benford", "Posterior")
    , col = c(grey(0.7), "black")
    , pch = c(16, 21)
    , pt.bg = c(NULL, "white")
    , lty = c("22", "solid")
    , lwd = c(1.25, 1)
    , bty = "n"
    , pt.cex = c(0.8, 1.5)
    , title = "Distribution"
    , seg.len = 1.5
  )