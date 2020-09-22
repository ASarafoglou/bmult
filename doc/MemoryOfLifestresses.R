## -----------------------------------------------------------------------------
library(multibridge)
data(lifestresses)
lifestresses

# visualize the data
op <- par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
  plot(lifestresses$month, lifestresses$stress.freq, col = "black", 
       pch = 1, cex = 2, type="b",
       xlim = c(0, 18), ylim = c(0, 20), 
       ylab = "", xlab = "", axes = FALSE)
  axis(1, at = seq(0, 18, length.out = 10), label = seq(0, 18, length.out = 10))
  axis(2, at = seq(0, 20, length.out = 5), label = seq(0, 20, length.out = 5), las = 1) 
  par(las = 0)
  mtext("Months Before Interview", side = 1, line = 2.5, cex = 1.5)
  mtext("Participants Reporting a\nNegative Life Event", side = 2, line = 3.0, cex = 1.5)

## -----------------------------------------------------------------------------
categories <- lifestresses$month
x          <- lifestresses$stress.freq

# Prior specification 
# We assign a uniform Dirichlet distribution, that is, we set all concentration parameters to 1
alpha <- rep(1, 18)

# Test the following restricted Hypothesis:
# Hr: month1 > month2 > ... > month18 
Hr   <- paste0(1:18, collapse=">"); Hr

## ---- cache = TRUE------------------------------------------------------------
ineq_results <- multBayesInformed(categories, Hr=Hr, a=alpha, counts=x, bf_type = 'BFre', seed = 2020)

## -----------------------------------------------------------------------------
summary(ineq_results)

## -----------------------------------------------------------------------------
ineq_results$bf_list
ineq_bayesfactors <- ineq_results$bf_list$bf

## ---- cache = TRUE------------------------------------------------------------
eq_results      <- multibridge::multBfEquality(a=alpha, counts=x)
eq_bayesfactors <- eq_results$bf

BFr0 <- ineq_bayesfactors[['BFre']] * eq_bayesfactors[['BFe0']]; BFr0

