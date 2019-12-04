#include <Rcpp.h>
#include "vwmvp.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cint_fun_J_RNplus(
    NumericVector x,
    NumericVector pars, 
    double radian) {
  
  int i;
  int lx = x.size();
  double kappa;
  double kc;
  NumericVector out(lx);
  
  for (i = 0; i < lx; i++) {
    
    if (i > 1e4) {
      kappa = x[i];
    } else {
    kappa = cKappaFromJ(x[i]);
    }
    
    kc = sqrt(pow(pars[2], 2) + pow(kappa, 2) + 2 * pars[2]*kappa*cos(radian));
    // Rcout << pow(pars[3], 2) << std::endl;
    // Rcout << kc << std::endl;
    
    out[i] = R::dgamma(x[i], pars[0]/pars[1], pars[1], 0);
    
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

