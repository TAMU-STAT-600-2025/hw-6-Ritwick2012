
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

library(testthat)
set.seed(123)

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################

test_that("soft thresholding C++ vs R", {
  # Input 1
  a1 <- 2.5
  lambda1 <- 1
  expect_equal(soft(a1, lambda1), soft_c(a1, lambda1))
  
  
  # Input 2
  a2 <- -0.8
  lambda2 <- 0.5
  expect_equal(soft(a2, lambda2), soft_c(a2, lambda2))
})

# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################

test_that("lasso objective C++ vs R", {
  # Generate data
  X <- matrix(rnorm(10 * 3), 10, 3)
  Y <- rnorm(10)
  
  # Input 1
  beta1 <- c(0.5, -0.2, 0)
  lambda1 <- 0.3
  expect_equal(lasso(X, Y, beta1, lambda1), lasso_c(X, Y, beta1, lambda1))
  
  
  # Input 2
  beta2 <- c(0.9, 0.1, 4)
  lambda2 <- 0.8
  expect_equal(lasso(X, Y, beta2, lambda2), lasso_c(X, Y, beta2, lambda2))
})

# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################

#Example 1 Simple case (n > p, small lambda)
n <- 20; p <- 5
X <- matrix(rnorm(n * p), n, p)
Y <- rnorm(n)

std <- standardizeXY(X, Y)
Xtilde <- std$Xtilde
Ytilde <- std$Ytilde
lambda <- 0.1

beta_r <- fitLASSOstandardized(Xtilde, Ytilde, lambda)$beta
beta_c <- fitLASSOstandardized_c(Xtilde, Ytilde, lambda, rep(0, p))

test_that("C++ vs R", {
  expect_equal(beta_r,as.vector(beta_c))
})


#Example 2: Lambda = 0 
n <- 25; p <- 8
X <- matrix(rnorm(n * p), n, p)
Y <- rnorm(n)

std <- standardizeXY(X, Y)
Xtilde <- std$Xtilde
Ytilde <- std$Ytilde
lambda <- 0

beta_r <- fitLASSOstandardized(Xtilde, Ytilde, lambda)$beta
beta_c <- fitLASSOstandardized_c(Xtilde, Ytilde, lambda, rep(0, p))

test_that("OLS case (lambda = 0): C++ vs R", {
  expect_equal(beta_r, as.vector(beta_c))
})


#Example 3: High-dimensional case (p>n)
n <- 10; p <- 15
X <- matrix(rnorm(n * p), n, p)
Y <- rnorm(n)

std <- standardizeXY(X, Y)
Xtilde <- std$Xtilde
Ytilde <- std$Ytilde
lambda <- 0.05

beta_r <- fitLASSOstandardized(Xtilde, Ytilde, lambda)$beta
beta_c <- fitLASSOstandardized_c(Xtilde, Ytilde, lambda, rep(0, p))

test_that("High-dimensional case: C++ vs R", {
  expect_equal(beta_r, as.vector(beta_c))
})

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

library(microbenchmark)
set.seed(20)
n <- 100
p <- 20
X <- matrix(rnorm(n*p), n, p)
Y <- rnorm(n)

# Standardize
std <- standardizeXY(X, Y)
Xtilde <- std$Xtilde
Ytilde <- std$Ytilde
beta_start <- rep(0, p)
lambda <- 0.1

# Run microbenchmark
mbm <- microbenchmark(
  R_version = fitLASSOstandardized(Xtilde, Ytilde, lambda, beta_start, eps = 1e-6),
  Cpp_version = fitLASSOstandardized_c(Xtilde, Ytilde, lambda, beta_start, eps = 1e-6),
  times = 50L
)

print(mbm)

# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Example 1 Simple case (n > p)
n <- 20; p <- 5
X <- matrix(rnorm(n * p), n, p)
Y <- rnorm(n)

std <- standardizeXY(X, Y)
Xtilde <- std$Xtilde
Ytilde <- std$Ytilde

lambda_seq <- c(0.5, 0.1, 0.05)

beta_r <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq)$beta_mat
beta_c <- fitLASSOstandardized_seq_c(Xtilde, Ytilde, lambda_seq)

test_that("n > p: C++ vs R sequence", {
  expect_equal(beta_r, beta_c)
})


# Example 2: High-dimensional p > n
n <- 10; p <- 15
X <- matrix(rnorm(n * p), n, p)
Y <- rnorm(n)

std <- standardizeXY(X, Y)
Xtilde <- std$Xtilde
Ytilde <- std$Ytilde

lambda_seq <- c(0.2, 0.1, 0.05)

beta_r <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq)$beta_mat
beta_c <- fitLASSOstandardized_seq_c(Xtilde, Ytilde, lambda_seq)

test_that("High-dimensional p > n: C++ vs R sequence", {
  expect_equal(beta_r, beta_c, tolerance = 1e-3)
})


# Example 3: Highly correlated predictors
n <- 30; p <- 10
base <- rnorm(n)
X <- matrix(0, n, p)
for (j in 1:p) X[, j] <- base + rnorm(n, sd = 0.01)  # highly correlated
Y <- rnorm(n)
std <- standardizeXY(X, Y)
Xtilde <- std$Xtilde
Ytilde <- std$Ytilde
lambda_seq <- c(0.2, 0.1)

beta_r3 <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq)$beta_mat
beta_c3 <- fitLASSOstandardized_seq_c(Xtilde, Ytilde, lambda_seq)

test_that("fitLASSOstandardized_seq: highly correlated predictors", {
  expect_equal(beta_r3, beta_c3)
})

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

# Generate example data
n <- 100
p <- 20
X <- matrix(rnorm(n*p), n, p)
Y <- rnorm(n)


# Standardize
std <- standardizeXY(X, Y)
Xtilde <- std$Xtilde
Ytilde <- std$Ytilde

# Define lambda sequence
lambda_seq <- seq(0.5, 0.01, length.out = 10)

# Run microbenchmark
mbm_seq <- microbenchmark(
  R_version_seq = fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq, eps = 1e-6),
  Cpp_version_seq = fitLASSOstandardized_seq_c(Xtilde, Ytilde, lambda_seq, eps = 1e-6),
  times = 50L
)

print(mbm_seq)

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)


