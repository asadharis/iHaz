#include <Rcpp.h>

using namespace Rcpp;

inline double do_myf( const NumericVector& x ) {
return sum(x);
 }
                                 
                                 // [[Rcpp::export]]
                                 List myf( List& X ) {
                                 
                                 return lapply( X, do_myf );
                                 
                                 }