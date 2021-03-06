// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fitGrove
Rcpp::List fitGrove(arma::mat D, arma::mat X, arma::vec p, arma::vec tau_par, arma::vec eta_par, arma::vec gamma_par, arma::vec init_state, double nu, double sigma0, double alpha, double beta, int n_samp, int transition_mode);
RcppExport SEXP grove_fitGrove(SEXP DSEXP, SEXP XSEXP, SEXP pSEXP, SEXP tau_parSEXP, SEXP eta_parSEXP, SEXP gamma_parSEXP, SEXP init_stateSEXP, SEXP nuSEXP, SEXP sigma0SEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP n_sampSEXP, SEXP transition_modeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau_par(tau_parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_par(eta_parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_par(gamma_parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init_state(init_stateSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_samp(n_sampSEXP);
    Rcpp::traits::input_parameter< int >::type transition_mode(transition_modeSEXP);
    __result = Rcpp::wrap(fitGrove(D, X, p, tau_par, eta_par, gamma_par, init_state, nu, sigma0, alpha, beta, n_samp, transition_mode));
    return __result;
END_RCPP
}
// fitGroveML
Rcpp::List fitGroveML(arma::mat D, arma::mat X, arma::vec p, arma::vec tau_par, arma::vec eta_par, arma::vec gamma_par, arma::vec init_state, double nu, double sigma0, double alpha, double beta, int transition_mode);
RcppExport SEXP grove_fitGroveML(SEXP DSEXP, SEXP XSEXP, SEXP pSEXP, SEXP tau_parSEXP, SEXP eta_parSEXP, SEXP gamma_parSEXP, SEXP init_stateSEXP, SEXP nuSEXP, SEXP sigma0SEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP transition_modeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau_par(tau_parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_par(eta_parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_par(gamma_parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init_state(init_stateSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type transition_mode(transition_modeSEXP);
    __result = Rcpp::wrap(fitGroveML(D, X, p, tau_par, eta_par, gamma_par, init_state, nu, sigma0, alpha, beta, transition_mode));
    return __result;
END_RCPP
}
