// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;


//This function builds the large D and D_star matrices
//which contain the terms d_{i,j} and d_{i,j}*
// [[Rcpp::export]]
List cpp_buildMatrixD(NumericVector ts, NumericVector as, NumericVector bs, NumericVector zs){
  int n = ts.size(), k = zs.size();

  //arma::colvec t(ts.begin(), n, false);
  //arma::colvec a(as.begin(), n, false);
  //arma::colvec b(bs.begin(), n, false);
  //arma::colvec z(zs.begin(), k, false);

  arma::mat D = arma::zeros(n,k-1);
  arma::mat Dstar = arma::zeros(n,k-1);


  for(int i = 0; i < n; i = i+1 )
  {
    for(int j = 1; j< k; j = j+1)
    {
      if(zs(j-1)>= as(i) && zs(j)<= bs(i))
      {
        Dstar(i,j-1) = 1;
      }
      if(zs(j-1)>= ts(i) && zs(j)<= as(i))
      {
        D(i,j-1) = 1;
      }
    }
  }

  return List::create(Named("D") = D,
  Named("Dstar") = Dstar);
}

// [[Rcpp::export]]
double cpp_logLikeRestricted(NumericVector lambda, NumericMatrix dsZ, NumericMatrix dZ){
 arma::colvec lam(lambda.begin(), lambda.size(), false);
 int r = dsZ.nrow(), c = dsZ.ncol();
 arma::mat dsz(dsZ.begin(), r, c, false);
 arma::mat dz(dZ.begin(), r, c, false);

 double result = accu(log(1-exp(-trans(lam)*dsz)) - trans(lam)*dz);

 return result;
}
// [[Rcpp::export]]
NumericVector cpp_nablaLogLikeRestrict(NumericVector lambda, NumericMatrix dsZ, NumericMatrix dZ){
 arma::colvec lam(lambda.begin(), lambda.size(), false);

 int r = dsZ.nrow(), c = dsZ.ncol();
 arma::mat dsz(dsZ.begin(), r, c, false);
 arma::mat dz(dZ.begin(), r, c, false);

 arma::mat res =  1/(exp(trans(lam)*dsz)-1);
 res = dsz*trans(res) - sum(dz, 1);
 return NumericVector(res.begin(),res.end());
}


// [[Rcpp::export]]
NumericVector cpp_newLambda(NumericVector pI, NumericVector lambda,NumericMatrix dsZ,NumericMatrix dZ){
  arma::colvec lam(lambda.begin(), lambda.size(), false);
  arma::colvec pi(pI.begin(), pI.size(), false);

  int r = dsZ.nrow(), c = dsZ.ncol();
  arma::mat dsz(dsZ.begin(), r, c, false);
  arma::mat dz(dZ.begin(), r, c, false);

  double alpha = 2.0;
  bool foundAlpha = false;

  while(!foundAlpha)
  {
    alpha =  alpha/2;
    NumericVector temp = cpp_nablaLogLikeRestrict(lambda,dsZ, dZ);
    arma::rowvec temp2(temp.begin(), temp.size(), false);

    foundAlpha = cpp_logLikeRestricted(lambda+alpha*pI, dsZ, dZ) -
        cpp_logLikeRestricted(lambda, dsZ, dZ)>=
        (alpha/4)*as_scalar(temp2*pi);
  }

  arma::colvec res = lam+alpha*pi;

  return NumericVector(res.begin(),res.end());
}

// [[Rcpp::export]]
NumericVector cpp_solve(NumericVector ys, NumericMatrix xs){
  int r = xs.nrow(), c = xs.ncol();

  arma::colvec y(ys.begin(), ys.size(), false);
  arma::mat x(xs.begin(), r,c,false);

  arma::colvec res = arma::solve(x,y);
  return NumericVector(res.begin(),res.end());
}

// [[Rcpp::export]]
NumericVector cpp_findLocalMax(NumericVector x, double tol){
  NumericVector index;
  int n  = x.size();

  for(int i = 1; i< n-1; i = i+1){
    if(x(i)>= x(i-1) && x(i)>= x(i+1) && x(i)> tol){
      index.push_back(i+1);
    }
  }
  return index;
}

// [[Rcpp::export]]
arma::mat cpp_scale(arma::mat matrix, arma::rowvec vec){
  matrix.each_row() %= vec;
  return matrix;
}


//THESE FUNCTONS WERE OBTAINED FROM THE EXAMPLES FILE FOR IMPLEMENTING A FAST TAPPLY
template <int RTYPE>
IntegerVector fast_factor_template( const Vector<RTYPE>& x ) {
    Vector<RTYPE> levs = sort_unique(x);
    IntegerVector out = match(x, levs);
    out.attr("levels") = as<CharacterVector>(levs);
    out.attr("class") = "factor";
    return out;
}

// [[Rcpp::export]]
SEXP fast_factor( SEXP x ) {
    switch( TYPEOF(x) ) {
    case INTSXP: return fast_factor_template<INTSXP>(x);
    case REALSXP: return fast_factor_template<REALSXP>(x);
    case STRSXP: return fast_factor_template<STRSXP>(x);
    }
    return R_NilValue;
}
