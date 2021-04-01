context("bridge sampling: get correct BF estimates")

input_samples <- cbind(c(0.740, 0.742, 0.742, 0.742, 0.741, 0.744, 0.744, 0.743, 0.741, 0.741,
                         0.748, 0.749, 0.750, 0.751, 0.752, 0.748, 0.745, 0.743, 0.741, 0.742),
                       c(0.260, 0.258, 0.258, 0.258, 0.259, 0.256, 0.256, 0.257, 0.259, 0.259, 
                         0.252, 0.251, 0.250, 0.249, 0.248, 0.252, 0.255, 0.257, 0.259, 0.258))
Hr            <- c('1  > 2')
counts        <- c(416, 140)
n             <- c(500, 200)
bridge_output <- mult_bf_inequality(samples=input_samples, Hr=Hr, x = counts, seed = 2)

test_that("mult: yields equal logml estimate with equal samples", {

  expect_equal(bridge_output$logml, -1.13495, tolerance = 0.002)
  
})

test_that("mult: yields equal RE2", {
  
  expect_equal(bridge_output$error_measures$re2, 0.793, tolerance = 0.002)
  
})

test_that("mult: prints equal hypothesis", {
  
  expect_equal(bridge_output$hyp, c("theta1", ">", "theta2"))
  
})


bridge_output <- binom_bf_inequality(samples=input_samples, Hr=Hr, x = counts, n=n, seed = 2)
test_that("binom: yields equal logml estimate with equal samples", {
  
  expect_equal(bridge_output$logml, -106.6986, tolerance = 0.002)
  
})

test_that("binom: yields equal RE2", {
  
  expect_equal(bridge_output$error_measures$re2, 0.8710915, tolerance = 0.002)
  
})

test_that("binom: prints equal hypothesis", {
  
  expect_equal(bridge_output$hyp, c("theta1", ">", "theta2"))
  
})
