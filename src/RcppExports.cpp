// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// computeCoxDirection
arma::mat computeCoxDirection(arma::vec& x, int comp, int n_points, int verbose);
RcppExport SEXP _opdesmixr_computeCoxDirection(SEXP xSEXP, SEXP compSEXP, SEXP n_pointsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type comp(compSEXP);
    Rcpp::traits::input_parameter< int >::type n_points(n_pointsSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(computeCoxDirection(x, comp, n_points, verbose));
    return rcpp_result_gen;
END_RCPP
}
// mixtureCoordinateExchangeGaussian
Rcpp::List mixtureCoordinateExchangeGaussian(arma::mat X_orig, int order, int n_cox_points, int max_it, int verbose);
RcppExport SEXP _opdesmixr_mixtureCoordinateExchangeGaussian(SEXP X_origSEXP, SEXP orderSEXP, SEXP n_cox_pointsSEXP, SEXP max_itSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_orig(X_origSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< int >::type n_cox_points(n_cox_pointsSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(mixtureCoordinateExchangeGaussian(X_orig, order, n_cox_points, max_it, verbose));
    return rcpp_result_gen;
END_RCPP
}
// getXsMNL
arma::mat getXsMNL(arma::cube& X, int s);
RcppExport SEXP _opdesmixr_getXsMNL(SEXP XSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(getXsMNL(X, s));
    return rcpp_result_gen;
END_RCPP
}
// getUsMNL
arma::vec getUsMNL(arma::cube& X, arma::vec& beta, int s, arma::mat& Xs);
RcppExport SEXP _opdesmixr_getUsMNL(SEXP XSEXP, SEXP betaSEXP, SEXP sSEXP, SEXP XsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Xs(XsSEXP);
    rcpp_result_gen = Rcpp::wrap(getUsMNL(X, beta, s, Xs));
    return rcpp_result_gen;
END_RCPP
}
// getPsMNL
arma::vec getPsMNL(arma::cube& X, arma::vec& beta, int s, arma::mat& Xs);
RcppExport SEXP _opdesmixr_getPsMNL(SEXP XSEXP, SEXP betaSEXP, SEXP sSEXP, SEXP XsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Xs(XsSEXP);
    rcpp_result_gen = Rcpp::wrap(getPsMNL(X, beta, s, Xs));
    return rcpp_result_gen;
END_RCPP
}
// getInformationMatrixMNL
arma::mat getInformationMatrixMNL(arma::cube& X, arma::vec& beta);
RcppExport SEXP _opdesmixr_getInformationMatrixMNL(SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(getInformationMatrixMNL(X, beta));
    return rcpp_result_gen;
END_RCPP
}
// getLogDEfficiencyMNL
double getLogDEfficiencyMNL(arma::cube& X, arma::vec& beta, int verbose);
RcppExport SEXP _opdesmixr_getLogDEfficiencyMNL(SEXP XSEXP, SEXP betaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(getLogDEfficiencyMNL(X, beta, verbose));
    return rcpp_result_gen;
END_RCPP
}
// findBestCoxDirMNL
arma::cube findBestCoxDirMNL(arma::mat& cox_dir, arma::cube& X_in, arma::vec& beta, int k, int s, double log_d_eff_best, int verbose);
RcppExport SEXP _opdesmixr_findBestCoxDirMNL(SEXP cox_dirSEXP, SEXP X_inSEXP, SEXP betaSEXP, SEXP kSEXP, SEXP sSEXP, SEXP log_d_eff_bestSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type cox_dir(cox_dirSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type X_in(X_inSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type log_d_eff_best(log_d_eff_bestSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(findBestCoxDirMNL(cox_dir, X_in, beta, k, s, log_d_eff_best, verbose));
    return rcpp_result_gen;
END_RCPP
}
// mixtureCoordinateExchangeMNL
Rcpp::List mixtureCoordinateExchangeMNL(arma::cube X_orig, arma::vec beta, int n_cox_points, int max_it, int verbose);
RcppExport SEXP _opdesmixr_mixtureCoordinateExchangeMNL(SEXP X_origSEXP, SEXP betaSEXP, SEXP n_cox_pointsSEXP, SEXP max_itSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type X_orig(X_origSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_cox_points(n_cox_pointsSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(mixtureCoordinateExchangeMNL(X_orig, beta, n_cox_points, max_it, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_opdesmixr_computeCoxDirection", (DL_FUNC) &_opdesmixr_computeCoxDirection, 4},
    {"_opdesmixr_mixtureCoordinateExchangeGaussian", (DL_FUNC) &_opdesmixr_mixtureCoordinateExchangeGaussian, 5},
    {"_opdesmixr_getXsMNL", (DL_FUNC) &_opdesmixr_getXsMNL, 2},
    {"_opdesmixr_getUsMNL", (DL_FUNC) &_opdesmixr_getUsMNL, 4},
    {"_opdesmixr_getPsMNL", (DL_FUNC) &_opdesmixr_getPsMNL, 4},
    {"_opdesmixr_getInformationMatrixMNL", (DL_FUNC) &_opdesmixr_getInformationMatrixMNL, 2},
    {"_opdesmixr_getLogDEfficiencyMNL", (DL_FUNC) &_opdesmixr_getLogDEfficiencyMNL, 3},
    {"_opdesmixr_findBestCoxDirMNL", (DL_FUNC) &_opdesmixr_findBestCoxDirMNL, 7},
    {"_opdesmixr_mixtureCoordinateExchangeMNL", (DL_FUNC) &_opdesmixr_mixtureCoordinateExchangeMNL, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_opdesmixr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}