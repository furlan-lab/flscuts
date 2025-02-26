
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
  const int nv = j.size();
  const int nm = rm.size();
  Rcpp::NumericVector rv(nm);
  Rcpp::NumericVector rit(nm);
  int current;
  // Calculate RowVars Initial
  for (int i = 0; i < nv; ++i) {
    current = j(i) - 1;
    rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
    rit(current) = rit(current) + 1;
  }
  // Calculate Remainder Variance
  for (int i = 0; i < nm; ++i) {
    rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
  }
  rv = rv / (n - 1);
  return(rv);
}
