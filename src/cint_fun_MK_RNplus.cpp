#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector cint_fun_MK_RNplus(
    NumericVector x,
    NumericVector pars, 
    double radian) {
  
  int i;
  int lx = x.size();
  double kappa;
  double kc;
  NumericVector out(lx);
  
  for (i = 0; i < lx; i++) {
    
    kappa = x[i];
    
    kc = sqrt(pow(pars[2], 2) + pow(kappa, 2) + 2 * pars[2]*kappa*cos(radian));
    // Rcout << pow(pars[3], 2) << std::endl;
    // Rcout << kc << std::endl;
    
    out[i] = R::dgamma(kappa, pars[0]/pars[1], pars[1], 0);
      
      if (out[i] != 0) {
        out[i] = out[i] * 
      ((R::bessel_i(kc,0,2) / 
      (2*PI*R::bessel_i(kappa,0,2) * 
      R::bessel_i(pars[2],0,2))) *
      exp(kc - (kappa + pars[2])));
      }
    
  }
  return(out);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
