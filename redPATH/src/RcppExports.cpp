// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// calc_vec
Rcpp::NumericVector calc_vec(Rcpp::NumericVector v1, Rcpp::NumericVector v2);
RcppExport SEXP _redPATH_calc_vec(SEXP v1SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type v2(v2SEXP);
    rcpp_result_gen = Rcpp::wrap(calc_vec(v1, v2));
    return rcpp_result_gen;
END_RCPP
}
// KL
Rcpp::NumericMatrix KL(const Rcpp::NumericMatrix& x, const Rcpp::NumericMatrix& logx);
RcppExport SEXP _redPATH_KL(SEXP xSEXP, SEXP logxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type logx(logxSEXP);
    rcpp_result_gen = Rcpp::wrap(KL(x, logx));
    return rcpp_result_gen;
END_RCPP
}
// ED2
Rcpp::NumericMatrix ED2(const Rcpp::NumericMatrix& x);
RcppExport SEXP _redPATH_ED2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(ED2(x));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP path_insertion_cost(SEXP, SEXP, SEXP);
RcppExport SEXP path_length_c(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_redPATH_calc_vec", (DL_FUNC) &_redPATH_calc_vec, 2},
    {"_redPATH_KL", (DL_FUNC) &_redPATH_KL, 2},
    {"_redPATH_ED2", (DL_FUNC) &_redPATH_ED2, 1},
    {"path_insertion_cost", (DL_FUNC) &path_insertion_cost, 3},
    {"path_length_c",       (DL_FUNC) &path_length_c,       2},
    {NULL, NULL, 0}
};

RcppExport void R_init_redPATH(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
