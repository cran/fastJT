// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// fastJT
Rcpp::List fastJT(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y, bool outTopNFlag, int outItemNo, bool standardized);
RcppExport SEXP fastJT_fastJT(SEXP XSEXP, SEXP YSEXP, SEXP outTopNFlagSEXP, SEXP outItemNoSEXP, SEXP standardizedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type outTopNFlag(outTopNFlagSEXP);
    Rcpp::traits::input_parameter< int >::type outItemNo(outItemNoSEXP);
    Rcpp::traits::input_parameter< bool >::type standardized(standardizedSEXP);
    __result = Rcpp::wrap(fastJT(X, Y, outTopNFlag, outItemNo, standardized));
    return __result;
END_RCPP
}
// fastJTmp
Rcpp::List fastJTmp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y, bool outTopNFlag, int numThreads, int outItemNo, bool standardized);
RcppExport SEXP fastJT_fastJTmp(SEXP XSEXP, SEXP YSEXP, SEXP outTopNFlagSEXP, SEXP numThreadsSEXP, SEXP outItemNoSEXP, SEXP standardizedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type outTopNFlag(outTopNFlagSEXP);
    Rcpp::traits::input_parameter< int >::type numThreads(numThreadsSEXP);
    Rcpp::traits::input_parameter< int >::type outItemNo(outItemNoSEXP);
    Rcpp::traits::input_parameter< bool >::type standardized(standardizedSEXP);
    __result = Rcpp::wrap(fastJTmp(X, Y, outTopNFlag, numThreads, outItemNo, standardized));
    return __result;
END_RCPP
}