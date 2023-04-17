// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

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
// getScheffeGaussian
arma::mat getScheffeGaussian(arma::mat& X, int order, int n_pv);
RcppExport SEXP _opdesmixr_getScheffeGaussian(SEXP XSEXP, SEXP orderSEXP, SEXP n_pvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    rcpp_result_gen = Rcpp::wrap(getScheffeGaussian(X, order, n_pv));
    return rcpp_result_gen;
END_RCPP
}
// getOptCritValueGaussian
double getOptCritValueGaussian(arma::mat& X, int order, int q, int opt_crit, arma::mat& W, int n_pv);
RcppExport SEXP _opdesmixr_getOptCritValueGaussian(SEXP XSEXP, SEXP orderSEXP, SEXP qSEXP, SEXP opt_critSEXP, SEXP WSEXP, SEXP n_pvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type opt_crit(opt_critSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    rcpp_result_gen = Rcpp::wrap(getOptCritValueGaussian(X, order, q, opt_crit, W, n_pv));
    return rcpp_result_gen;
END_RCPP
}
// changeIngredientDesignCoxGaussian
void changeIngredientDesignCoxGaussian(double theta, arma::mat& X, int i, int j, int n_pv);
RcppExport SEXP _opdesmixr_changeIngredientDesignCoxGaussian(SEXP thetaSEXP, SEXP XSEXP, SEXP iSEXP, SEXP jSEXP, SEXP n_pvSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    changeIngredientDesignCoxGaussian(theta, X, i, j, n_pv);
    return R_NilValue;
END_RCPP
}
// efficiencyCoxScheffeGaussian
double efficiencyCoxScheffeGaussian(double theta, arma::mat& X, int i, int j, int order, int opt_crit, arma::mat& W, int n_pv);
RcppExport SEXP _opdesmixr_efficiencyCoxScheffeGaussian(SEXP thetaSEXP, SEXP XSEXP, SEXP iSEXP, SEXP jSEXP, SEXP orderSEXP, SEXP opt_critSEXP, SEXP WSEXP, SEXP n_pvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< int >::type opt_crit(opt_critSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    rcpp_result_gen = Rcpp::wrap(efficiencyCoxScheffeGaussian(theta, X, i, j, order, opt_crit, W, n_pv));
    return rcpp_result_gen;
END_RCPP
}
// mixtureCoordinateExchangeGaussian
Rcpp::List mixtureCoordinateExchangeGaussian(const arma::mat X_orig, int order, int max_it, int verbose, int opt_crit, arma::mat W, int opt_method, double lower, double upper, double tol, int n_cox_points, int n_pv);
RcppExport SEXP _opdesmixr_mixtureCoordinateExchangeGaussian(SEXP X_origSEXP, SEXP orderSEXP, SEXP max_itSEXP, SEXP verboseSEXP, SEXP opt_critSEXP, SEXP WSEXP, SEXP opt_methodSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP tolSEXP, SEXP n_cox_pointsSEXP, SEXP n_pvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X_orig(X_origSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type opt_crit(opt_critSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type opt_method(opt_methodSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type n_cox_points(n_cox_pointsSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    rcpp_result_gen = Rcpp::wrap(mixtureCoordinateExchangeGaussian(X_orig, order, max_it, verbose, opt_crit, W, opt_method, lower, upper, tol, n_cox_points, n_pv));
    return rcpp_result_gen;
END_RCPP
}
// getXsMNL
arma::mat getXsMNL(arma::cube& X, int s, int order, int n_pv, bool no_choice);
RcppExport SEXP _opdesmixr_getXsMNL(SEXP XSEXP, SEXP sSEXP, SEXP orderSEXP, SEXP n_pvSEXP, SEXP no_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    Rcpp::traits::input_parameter< bool >::type no_choice(no_choiceSEXP);
    rcpp_result_gen = Rcpp::wrap(getXsMNL(X, s, order, n_pv, no_choice));
    return rcpp_result_gen;
END_RCPP
}
// getUsMNL
arma::vec getUsMNL(arma::cube& X, arma::vec& beta, int s, arma::mat& Xs, bool transform_beta, int n_pv);
RcppExport SEXP _opdesmixr_getUsMNL(SEXP XSEXP, SEXP betaSEXP, SEXP sSEXP, SEXP XsSEXP, SEXP transform_betaSEXP, SEXP n_pvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< bool >::type transform_beta(transform_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    rcpp_result_gen = Rcpp::wrap(getUsMNL(X, beta, s, Xs, transform_beta, n_pv));
    return rcpp_result_gen;
END_RCPP
}
// getPsMNL
arma::vec getPsMNL(arma::cube& X, arma::vec& beta, int s, arma::mat& Xs, bool transform_beta, int n_pv);
RcppExport SEXP _opdesmixr_getPsMNL(SEXP XSEXP, SEXP betaSEXP, SEXP sSEXP, SEXP XsSEXP, SEXP transform_betaSEXP, SEXP n_pvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< bool >::type transform_beta(transform_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    rcpp_result_gen = Rcpp::wrap(getPsMNL(X, beta, s, Xs, transform_beta, n_pv));
    return rcpp_result_gen;
END_RCPP
}
// getInformationMatrixMNL
arma::mat getInformationMatrixMNL(arma::cube& X, arma::vec& beta, int order, bool transform_beta, int n_pv, bool no_choice);
RcppExport SEXP _opdesmixr_getInformationMatrixMNL(SEXP XSEXP, SEXP betaSEXP, SEXP orderSEXP, SEXP transform_betaSEXP, SEXP n_pvSEXP, SEXP no_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type transform_beta(transform_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    Rcpp::traits::input_parameter< bool >::type no_choice(no_choiceSEXP);
    rcpp_result_gen = Rcpp::wrap(getInformationMatrixMNL(X, beta, order, transform_beta, n_pv, no_choice));
    return rcpp_result_gen;
END_RCPP
}
// getBayesianInformationMatrixMNL
arma::mat getBayesianInformationMatrixMNL(arma::cube& X, arma::mat& beta_mat, int order, bool transform_beta, int n_pv, bool no_choice);
RcppExport SEXP _opdesmixr_getBayesianInformationMatrixMNL(SEXP XSEXP, SEXP beta_matSEXP, SEXP orderSEXP, SEXP transform_betaSEXP, SEXP n_pvSEXP, SEXP no_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type transform_beta(transform_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    Rcpp::traits::input_parameter< bool >::type no_choice(no_choiceSEXP);
    rcpp_result_gen = Rcpp::wrap(getBayesianInformationMatrixMNL(X, beta_mat, order, transform_beta, n_pv, no_choice));
    return rcpp_result_gen;
END_RCPP
}
// getOptCritValueMNL
double getOptCritValueMNL(arma::cube& X, arma::mat& beta_mat, int verbose, int opt_crit, arma::mat& W, int order, bool transform_beta, int n_pv, bool no_choice);
RcppExport SEXP _opdesmixr_getOptCritValueMNL(SEXP XSEXP, SEXP beta_matSEXP, SEXP verboseSEXP, SEXP opt_critSEXP, SEXP WSEXP, SEXP orderSEXP, SEXP transform_betaSEXP, SEXP n_pvSEXP, SEXP no_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type opt_crit(opt_critSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type transform_beta(transform_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    Rcpp::traits::input_parameter< bool >::type no_choice(no_choiceSEXP);
    rcpp_result_gen = Rcpp::wrap(getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W, order, transform_beta, n_pv, no_choice));
    return rcpp_result_gen;
END_RCPP
}
// findBestCoxDirMNLDiscrete
void findBestCoxDirMNLDiscrete(arma::mat& cox_dir, arma::cube& X, arma::mat& beta_mat, int k, int s, double opt_crit_value_best, int verbose, int opt_crit, arma::mat& W, int order, bool transform_beta);
RcppExport SEXP _opdesmixr_findBestCoxDirMNLDiscrete(SEXP cox_dirSEXP, SEXP XSEXP, SEXP beta_matSEXP, SEXP kSEXP, SEXP sSEXP, SEXP opt_crit_value_bestSEXP, SEXP verboseSEXP, SEXP opt_critSEXP, SEXP WSEXP, SEXP orderSEXP, SEXP transform_betaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type cox_dir(cox_dirSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type opt_crit_value_best(opt_crit_value_bestSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type opt_crit(opt_critSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type transform_beta(transform_betaSEXP);
    findBestCoxDirMNLDiscrete(cox_dir, X, beta_mat, k, s, opt_crit_value_best, verbose, opt_crit, W, order, transform_beta);
    return R_NilValue;
END_RCPP
}
// changeIngredientDesignMNL
void changeIngredientDesignMNL(double theta, arma::cube& X, int i, int j, int s, int n_pv);
RcppExport SEXP _opdesmixr_changeIngredientDesignMNL(SEXP thetaSEXP, SEXP XSEXP, SEXP iSEXP, SEXP jSEXP, SEXP sSEXP, SEXP n_pvSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    changeIngredientDesignMNL(theta, X, i, j, s, n_pv);
    return R_NilValue;
END_RCPP
}
// efficiencyCoxScheffeMNL
double efficiencyCoxScheffeMNL(double theta, arma::cube& X, arma::mat& beta_mat, int i, int j, int s, int opt_crit, arma::mat& W, int order, bool transform_beta, int n_pv, bool no_choice);
RcppExport SEXP _opdesmixr_efficiencyCoxScheffeMNL(SEXP thetaSEXP, SEXP XSEXP, SEXP beta_matSEXP, SEXP iSEXP, SEXP jSEXP, SEXP sSEXP, SEXP opt_critSEXP, SEXP WSEXP, SEXP orderSEXP, SEXP transform_betaSEXP, SEXP n_pvSEXP, SEXP no_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type opt_crit(opt_critSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type transform_beta(transform_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    Rcpp::traits::input_parameter< bool >::type no_choice(no_choiceSEXP);
    rcpp_result_gen = Rcpp::wrap(efficiencyCoxScheffeMNL(theta, X, beta_mat, i, j, s, opt_crit, W, order, transform_beta, n_pv, no_choice));
    return rcpp_result_gen;
END_RCPP
}
// findBestCoxDirMNLBrent
void findBestCoxDirMNLBrent(arma::cube& X, arma::mat& beta_mat, int i, int j, int s, int opt_crit, int order, arma::mat& W, double lower, double upper, double tol, int verbose, bool transform_beta, int n_pv, bool no_choice);
RcppExport SEXP _opdesmixr_findBestCoxDirMNLBrent(SEXP XSEXP, SEXP beta_matSEXP, SEXP iSEXP, SEXP jSEXP, SEXP sSEXP, SEXP opt_critSEXP, SEXP orderSEXP, SEXP WSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP tolSEXP, SEXP verboseSEXP, SEXP transform_betaSEXP, SEXP n_pvSEXP, SEXP no_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type opt_crit(opt_critSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type transform_beta(transform_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    Rcpp::traits::input_parameter< bool >::type no_choice(no_choiceSEXP);
    findBestCoxDirMNLBrent(X, beta_mat, i, j, s, opt_crit, order, W, lower, upper, tol, verbose, transform_beta, n_pv, no_choice);
    return R_NilValue;
END_RCPP
}
// mixtureCoordinateExchangeMNL
Rcpp::List mixtureCoordinateExchangeMNL(arma::cube X_orig, arma::mat beta_mat, int order, int max_it, int verbose, int opt_crit, arma::mat W, int opt_method, double lower, double upper, double tol, int n_cox_points, bool transform_beta, int n_pv, bool no_choice);
RcppExport SEXP _opdesmixr_mixtureCoordinateExchangeMNL(SEXP X_origSEXP, SEXP beta_matSEXP, SEXP orderSEXP, SEXP max_itSEXP, SEXP verboseSEXP, SEXP opt_critSEXP, SEXP WSEXP, SEXP opt_methodSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP tolSEXP, SEXP n_cox_pointsSEXP, SEXP transform_betaSEXP, SEXP n_pvSEXP, SEXP no_choiceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type X_orig(X_origSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type opt_crit(opt_critSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type opt_method(opt_methodSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type n_cox_points(n_cox_pointsSEXP);
    Rcpp::traits::input_parameter< bool >::type transform_beta(transform_betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_pv(n_pvSEXP);
    Rcpp::traits::input_parameter< bool >::type no_choice(no_choiceSEXP);
    rcpp_result_gen = Rcpp::wrap(mixtureCoordinateExchangeMNL(X_orig, beta_mat, order, max_it, verbose, opt_crit, W, opt_method, lower, upper, tol, n_cox_points, transform_beta, n_pv, no_choice));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_opdesmixr_computeCoxDirection", (DL_FUNC) &_opdesmixr_computeCoxDirection, 4},
    {"_opdesmixr_getScheffeGaussian", (DL_FUNC) &_opdesmixr_getScheffeGaussian, 3},
    {"_opdesmixr_getOptCritValueGaussian", (DL_FUNC) &_opdesmixr_getOptCritValueGaussian, 6},
    {"_opdesmixr_changeIngredientDesignCoxGaussian", (DL_FUNC) &_opdesmixr_changeIngredientDesignCoxGaussian, 5},
    {"_opdesmixr_efficiencyCoxScheffeGaussian", (DL_FUNC) &_opdesmixr_efficiencyCoxScheffeGaussian, 8},
    {"_opdesmixr_mixtureCoordinateExchangeGaussian", (DL_FUNC) &_opdesmixr_mixtureCoordinateExchangeGaussian, 12},
    {"_opdesmixr_getXsMNL", (DL_FUNC) &_opdesmixr_getXsMNL, 5},
    {"_opdesmixr_getUsMNL", (DL_FUNC) &_opdesmixr_getUsMNL, 6},
    {"_opdesmixr_getPsMNL", (DL_FUNC) &_opdesmixr_getPsMNL, 6},
    {"_opdesmixr_getInformationMatrixMNL", (DL_FUNC) &_opdesmixr_getInformationMatrixMNL, 6},
    {"_opdesmixr_getBayesianInformationMatrixMNL", (DL_FUNC) &_opdesmixr_getBayesianInformationMatrixMNL, 6},
    {"_opdesmixr_getOptCritValueMNL", (DL_FUNC) &_opdesmixr_getOptCritValueMNL, 9},
    {"_opdesmixr_findBestCoxDirMNLDiscrete", (DL_FUNC) &_opdesmixr_findBestCoxDirMNLDiscrete, 11},
    {"_opdesmixr_changeIngredientDesignMNL", (DL_FUNC) &_opdesmixr_changeIngredientDesignMNL, 6},
    {"_opdesmixr_efficiencyCoxScheffeMNL", (DL_FUNC) &_opdesmixr_efficiencyCoxScheffeMNL, 12},
    {"_opdesmixr_findBestCoxDirMNLBrent", (DL_FUNC) &_opdesmixr_findBestCoxDirMNLBrent, 15},
    {"_opdesmixr_mixtureCoordinateExchangeMNL", (DL_FUNC) &_opdesmixr_mixtureCoordinateExchangeMNL, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_opdesmixr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
