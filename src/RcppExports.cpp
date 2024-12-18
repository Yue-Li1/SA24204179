// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mean_variance_model
List mean_variance_model(NumericVector expected_returns, NumericMatrix covariance_matrix, NumericVector weights);
RcppExport SEXP _SA24204179_mean_variance_model(SEXP expected_returnsSEXP, SEXP covariance_matrixSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type expected_returns(expected_returnsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covariance_matrix(covariance_matrixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(mean_variance_model(expected_returns, covariance_matrix, weights));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA24204179_mean_variance_model", (DL_FUNC) &_SA24204179_mean_variance_model, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA24204179(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
