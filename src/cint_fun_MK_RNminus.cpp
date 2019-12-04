#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cint_fun_MK_RNminus(
    NumericVector x,
    NumericVector pars, 
    double radian) {
  
  int i;
  int lx = x.size();
  NumericVector out(lx);
  
  for (i = 0; i < lx; i++) {
    
    out[i] = R::dgamma(x[i], pars[0]/pars[1], pars[1], 0);
    
    if (out[i] != 0) {
    out[i] = out [i] * 
      (1 / (2*PI*R::bessel_i(x[i],0,2)) * 
      pow(exp(cos(radian) - 1),x[i]));
    }
    
  }
  return(out);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


