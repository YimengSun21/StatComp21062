#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler to generate a chain with target joint density
//' @param n parameter in Bernoulli distribution
//' @param a parameter in Beta distribution
//' @param b parameter in Beta distribution
//' @return a two-dimensional chain of size 10000
//' @examples
//' \dontrun{
//' sampleC = gibbsC(50, 30, 40)
//' plot(sampleC, cex = 0.5)
//' }
//' @importFrom stats rbinom rbeta
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int n, double a, double b) {
  NumericMatrix mat (10000, 2);
  mat(0,0) = 0;
  mat(0,1) = 0;
  for (int i = 1; i < 10000; i++) {
    int xt = rbinom(1, n, mat(i-1,1))[0];
    double yt = rbeta(1, xt + a, n - xt + b)[0];
    mat(i,0) = xt;
    mat(i,1) = yt;
  };
  return mat;
}
