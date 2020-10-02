## -----------------------------------------------------------------------------
library(multibridge)
data(peas)
peas

## -----------------------------------------------------------------------------
categories <- peas$peas
x          <- peas$counts

# Prior specification 
# We assign a uniform Dirichlet distribution, that is, we set all concentration parameters to 1
alpha <- c(1, 1, 1, 1)

# Test the following restricted Hypothesis:
# Hr: month1 > month2 > ... > month18 
Hr   <- c('roundYellow > wrinkledYellow = roundGreen > wrinkledGreen')

## ---- cache = TRUE------------------------------------------------------------
results <- multibridge::multBfInformed(categories, Hr=Hr, a=alpha, counts=x, bf_type = 'BFre', seed = 2020)

## -----------------------------------------------------------------------------
summary(results)
# summary(ineq_results)$equalities$bf

