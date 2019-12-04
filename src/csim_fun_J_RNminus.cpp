// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>
#include "vwmvp.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector csim_fun_J_RNminus(
    double x,
    NumericVector pars, 
    NumericVector baseradians) {
  
  int i;
  int radian;
  NumericVector J;
  double kappa;
  int lr = baseradians.size();
  NumericVector out(lr);
  double radians = 0;
  
  J = Rcpp::rgamma(x,pars[0]/pars[1], pars[1]);
  
  for (radian = 0; radian < lr; radian++){
    
    for (i = 0; i < x; i++) {
      
      kappa = cKappaFromJ(J[i]);
      
      radians = radians + 
        1/(2*PI*exp(-kappa)*boost::math::cyl_bessel_i(0, kappa)) *
          pow(exp(cos(baseradians[radian]) - 1 ),kappa);
      
    }
    
    out[radian] = radians;
    radians = 0;
  }
  
  return(out/sum(out));
}
