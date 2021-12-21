#include <Rcpp.h>
using namespace Rcpp;
//' @title A Gibbs sampler using Rcpp
//' @description A specific Gibbs sampler for 2 dimension (x,y) with bivariate density using Rcpp
//' @param m numbers of sample want to generate 
//' @param a the a parameter in Beta(x+a,n-x+b)
//' @param b the b parameter in Beta(x+a,n-x+b)
//' @param n the n parameter in Binomial(n,y)
//' @return a matrix including two random sample with dimension of (\code{m},2)
//' @examples
//' \dontrun{
//'   dir_cpp <- '/data/cenmin/statistical computing/Rcpp/'
//'   sourceCpp(paste0(dir_cpp,"GBSampler.cpp"))
//'   set.seed(12)
//'   m = 10000
//'   a = 2
//'   b = 3
//'   n = 5
//'   t1 <- proc.time()
//'   xC <- gbsamplerC(m,a,b,n)
//'   t2 <- proc.time()
//'   t <- t2-t1
//'   paste0(t[3][[1]],'s') 
//' } 
//' @export
// [[Rcpp::export]]
NumericMatrix gbsamplerC(int m, int a, int b, int n) {
  NumericMatrix x(m,2);
  x(0,0) = 2.5;
  x(0,1) = 0.5;
  for(int i=1; i < m; i++){
    x(i,0) = R::rbinom(n,x(i-1,1));
    x(i,1) = R::rbeta(x(i-1,0)+a,n-x(i-1,0)+b);
  }
  return x;
}


