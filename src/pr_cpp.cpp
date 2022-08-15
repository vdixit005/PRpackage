#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List pr_cpp(NumericVector f0, NumericMatrix U, NumericMatrix X, Function d, int N, NumericVector w) {
  int n = X.rows();
  int t = U.rows();
  NumericVector favg = 0.0 * f0;
  double Lavg = 0.0;
  NumericVector Davg = 0.0 * f0;
  for(int j = 0; j < N; j++){
    NumericVector f = f0 + 0.0;
    double L = 0.0;
    NumericVector D (t, 1.0);
    NumericMatrix kernel(n,t);
    for(int i = 0; i < n; i++){
      NumericVector num (t,0.0);
      NumericVector num2 (t,0.0);
      for(int k = 0; k < t; k++){
        kernel(i , k) = as<double>(d(X( i , _ ), U( k , _ )));
        num(k) = kernel(i , k)*f(k);
        num2(k) = kernel(i , k)*D(k);
      }
      double den = Rcpp::sum(num2) / t;
      L = L + log(den);
      f = (1.0 - w(i))*f + w(i)*num / den;
      D = D*(1.0 + w(i)*(kernel(i, _) / den - 1.0));
    }
    favg = (j)*favg/(j+1) + f/(j+1);
    Lavg = (j)*Lavg/(j+1) + L/(j+1);
    Davg = (j)*Davg/(j+1) + D/(j+1);
  }
  List final;
  final("f") = favg;
  final("L") = Lavg;
  final("D") = Davg;
  return final;
}

