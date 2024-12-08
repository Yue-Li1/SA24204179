## -----------------------------------------------------------------------------
admm_algorithm <- function(A, B, C, f, g, rho = 1, max_iter = 100, tol = 1e-6) {
  x <- matrix(0, nrow = ncol(A), ncol = 1)
  z <- matrix(0, nrow = ncol(B), ncol = 1)
  y <- matrix(0, nrow = nrow(A), ncol = 1)
  lambda <- matrix(0, nrow = nrow(A), ncol = 1)
  for (iter in 1:max_iter) {
    x <- solve(t(A) %*% A + rho * diag(nrow(A)), 
               t(A) %*% (C - B %*% z) + rho * (y - lambda / rho))
    z <- solve(t(B) %*% B + rho * diag(nrow(B)), 
               t(B) %*% (C - A %*% x) + rho * (y - lambda / rho))
    y <- x
    lambda <- lambda + rho * (A %*% x + B %*% z - C)
    if (sum((A %*% x + B %*% z - C)^2) < tol) {
      cat("Converged at iteration", iter, "\n")
      break
    }
  }
  return(list(x = x, z = z, lambda = lambda))
}
A <- matrix(c(1, 2, 3, 4), 2, 2)
B <- diag(2)
C <- c(1, 1)
f <- function(x) sum(x^2)
g <- function(z) sum(abs(z))
result <- admm_algorithm(A, B, C, f, g)
print(result$x)

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
List mean_variance_model(NumericVector expected_returns, NumericMatrix covariance_matrix, NumericVector weights) {
  int n = expected_returns.size();
  double mean = sum(expected_returns * weights);
  double variance = 0.0;
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n; ++j) {
     variance += weights[i] * weights[j] * covariance_matrix(i, j);
    }
  }
  return List::create(
    Named("mean") = mean,
    Named("variance") = variance
  );
}')
expected_returns <- c(0.1, 0.12, 0.15)
covariance_matrix <- matrix(c(0.01, 0.001, 0.002, 
                              0.001, 0.02, 0.003, 
                              0.002, 0.003, 0.03), nrow = 3)
weights <- c(0.4, 0.4, 0.2)
result <- mean_variance_model(expected_returns, covariance_matrix, weights)
print(result)

