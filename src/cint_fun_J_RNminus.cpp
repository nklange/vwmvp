#include <Rcpp.h>
#include"vwmvp.hpp"
using namespace Rcpp;

//' Rcpp core model functions for variable precision models
//'
//' @param x A numeric vector of kappa/J values
//' @param pars A numeric vector of parameter values: RN- c(setsize precision,tau); RN+ c(setsize precision,tau,kappa_r)
//' @param radian A double for error in radian
//' @export


// [[Rcpp::export]]
NumericVector cint_fun_J_RNminus(
    NumericVector x,
    NumericVector pars, 
    double radian) {
  
  int i;
  int lx = x.size();
  double kappa;
  NumericVector out(lx);
  
  for (i = 0; i < lx; i++) {
    
    if (i > 1e4) {
      kappa = x[i];
    } else {
    kappa = cKappaFromJ(x[i]);
    }
    
    
    out[i] = R::dgamma(x[i], pars[0]/pars[1], pars[1], 0);
    
    if (out[i] != 0) {
      out[i] = out[i] *
      (1 / (2*PI*R::bessel_i(kappa,0,2)) * 
      pow(exp(cos(radian) - 1),kappa));
    }
    
  }
  return(out);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


