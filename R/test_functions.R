# Load the necessary packages
library(testthat)

# Load simulation functions
# source("sim_functions.R")

# Define test cases

# Test opponentprocess() with default parameters
test_that("opponentprocess() with default parameters runs without errors", {
  expect_no_error(opponentprocess())
})

# Test opponentprocess() with custom parameters and plot utility function
test_that("opponentprocess() with custom parameters and utility plot runs without errors", {
  expect_no_error(opponentprocess(plot_utility = TRUE, k_apk = 0.01, k_bpk = 0.01, k_apd = 1, k_bpd = 0.01, k_H = 1, gamma_a = 0.5, lambda_b = 1, gamma_b = 0.7, infuse = 1))
})

# Test bode_plot() with default parameters
test_that("bode_plot() with default parameters runs without errors", {
  expect_no_error(bode_plot())
})

# Test bode_plot() with custom parameters
test_that("bode_plot() with custom parameters runs without errors", {
  expect_no_error(bode_plot(gamma_a = 0.2, gamma_b = c(0.5, 0.7), lambda_b = 1, k_apk = 0.005, k_bpk = 0.004, freq_interval = 0.0002, multiply = 40, plot_frequencies = c(0.0002, 0.006)))
})

# Run the tests
test_results <- test_dir(".", reporter = "progress")

# Print the test results
print(test_results)
