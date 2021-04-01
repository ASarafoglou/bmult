library('multibridge')

# Load Data
  data(journals)

# Since percentages are rounded to two decimal values, we round the
# articles with an error to obtain integer values
  x <- round(journals$articles_with_NHST  * 
               (journals$perc_articles_with_errors/100))
# Total number of articles
  n <- journals$articles_with_NHST

# Prior specification
# We assign a uniform beta distribution to each binomial proportion
  a <- rep(1, 8)
  b <- rep(1, 8)

# Specifying the informed Hypothesis
  Hr <- c('JAP , PS , JCCP , PLOS , DP , FP , JEPG < JPSP')

# Category labels
  journal_names <- journals$journal

# Execute the analysis
  results_H0_Hr <- binom_bf_informed(x = x, n = n, Hr = Hr, a = a, b = b,
                                     factor_levels = journal_names,
                                     bf_type = 'LogBFr0', seed = 2020)

# Summarize Bayes Factors
  LogBFr0 <- summary(results_H0_Hr)$bf
  LogBFe0 <- results_H0_Hr$bf_list$bf0_table[['LogBFe0']]
  LogBFre <- -results_H0_Hr$bf_list$bfr_table[['LogBFer']]
  
  bayes_factor_table <- data.frame(
    BFType = c('LogBFe0', 'LogBFr0', 'LogBFre'), 
    BF = c(LogBFe0, LogBFr0, LogBFre)
    )
  
# Create Plot
  plot(summary(results_H0_Hr), xlab = "Journal")