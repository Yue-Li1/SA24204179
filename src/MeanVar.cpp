#include <Rcpp.h>
using namespace Rcpp;
//' @title Mean-Variance Model Calculation.
//' @description Calculate the mean and variance of the portfolio.
//' @param expected_returns A numeric vector of expected returns for each asset.
//' @param covariance_matrix A matrix representing the covariance between assets.
//' @param weights A numeric vector of weights for each asset in the portfolio.
//' @return A list containing the mean and variance of the portfolio.
//' @examples
//' \dontrun{
//' expected_returns <- c(0.1, 0.12, 0.15)
//' covariance_matrix <- matrix(c(0.01, 0.001, 0.002, 
//'                               0.001, 0.02, 0.003, 
//'                               0.002, 0.003, 0.03), nrow = 3)
//' weights <- c(0.4, 0.4, 0.2)
//' result <- mean_variance_model(expected_returns, covariance_matrix, weights)
//' print(result)
//' }
//' @export
// [[Rcpp::export]]
List mean_variance_model(NumericVector expected_returns, NumericMatrix covariance_matrix, NumericVector weights) {
  int n = expected_returns.size();
  // Calculate the expected return of the portfolio
  double mean = sum(expected_returns * weights);
  // Calculate the variance of the portfolio
  double variance = 0.0;
  // Loop over all pairs of assets (i, j)
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n; ++j) {
     variance += weights[i] * weights[j] * covariance_matrix(i, j);  // Weighted covariance
    }
  }
  // Return a list containing the mean and variance of the portfolio
  return List::create(
    Named("mean") = mean,  // Named to clarify what the values represent
    Named("variance") = variance
  );
}
 