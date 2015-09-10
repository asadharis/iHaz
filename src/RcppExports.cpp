// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cpp_buildMatrixD
List cpp_buildMatrixD(NumericVector ts, NumericVector as, NumericVector bs, NumericVector zs);
RcppExport SEXP iHaz_cpp_buildMatrixD(SEXP tsSEXP, SEXP asSEXP, SEXP bsSEXP, SEXP zsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type ts(tsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type as(asSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    __result = Rcpp::wrap(cpp_buildMatrixD(ts, as, bs, zs));
    return __result;
END_RCPP
}
// cpp_logLikeRestricted
double cpp_logLikeRestricted(NumericVector lambda, NumericMatrix dsZ, NumericMatrix dZ);
RcppExport SEXP iHaz_cpp_logLikeRestricted(SEXP lambdaSEXP, SEXP dsZSEXP, SEXP dZSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dsZ(dsZSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dZ(dZSEXP);
    __result = Rcpp::wrap(cpp_logLikeRestricted(lambda, dsZ, dZ));
    return __result;
END_RCPP
}
// cpp_nablaLogLikeRestrict
NumericVector cpp_nablaLogLikeRestrict(NumericVector lambda, NumericMatrix dsZ, NumericMatrix dZ);
RcppExport SEXP iHaz_cpp_nablaLogLikeRestrict(SEXP lambdaSEXP, SEXP dsZSEXP, SEXP dZSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dsZ(dsZSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dZ(dZSEXP);
    __result = Rcpp::wrap(cpp_nablaLogLikeRestrict(lambda, dsZ, dZ));
    return __result;
END_RCPP
}
// cpp_newLambda
NumericVector cpp_newLambda(NumericVector pI, NumericVector lambda, NumericMatrix dsZ, NumericMatrix dZ);
RcppExport SEXP iHaz_cpp_newLambda(SEXP pISEXP, SEXP lambdaSEXP, SEXP dsZSEXP, SEXP dZSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type pI(pISEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dsZ(dsZSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dZ(dZSEXP);
    __result = Rcpp::wrap(cpp_newLambda(pI, lambda, dsZ, dZ));
    return __result;
END_RCPP
}
// cpp_solve
NumericVector cpp_solve(NumericVector ys, NumericMatrix xs);
RcppExport SEXP iHaz_cpp_solve(SEXP ysSEXP, SEXP xsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xs(xsSEXP);
    __result = Rcpp::wrap(cpp_solve(ys, xs));
    return __result;
END_RCPP
}
// cpp_findLocalMax
NumericVector cpp_findLocalMax(NumericVector x, double tol);
RcppExport SEXP iHaz_cpp_findLocalMax(SEXP xSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    __result = Rcpp::wrap(cpp_findLocalMax(x, tol));
    return __result;
END_RCPP
}
// cpp_scale
arma::mat cpp_scale(arma::mat matrix, arma::rowvec vec);
RcppExport SEXP iHaz_cpp_scale(SEXP matrixSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type vec(vecSEXP);
    __result = Rcpp::wrap(cpp_scale(matrix, vec));
    return __result;
END_RCPP
}
// fast_factor
SEXP fast_factor(SEXP x);
RcppExport SEXP iHaz_fast_factor(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    __result = Rcpp::wrap(fast_factor(x));
    return __result;
END_RCPP
}
// myf
List myf(List& X);
RcppExport SEXP iHaz_myf(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List& >::type X(XSEXP);
    __result = Rcpp::wrap(myf(X));
    return __result;
END_RCPP
}
