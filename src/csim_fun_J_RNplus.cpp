// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>
#include "vwmvp.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector csim_fun_J_RNplus(
    double x,
    NumericVector pars, 
    NumericVector baseradians) {
  
  int i;
  int radian;
  double kappa;
  NumericVector J;
  double kc;
  int lr = baseradians.size();
  NumericVector out(lr);
  double radians = 0;
  
  J = Rcpp::rgamma(x,pars[0]/pars[1], pars[1]);
  
  for (radian = 0; radian < lr; radian++){
    
    for (i = 0; i < x; i++) {
      
      kappa = cKappaFromJ(J[i]);
      
      kc = sqrt(pow(pars[2], 2) + pow(kappa, 2) + 2 * pars[2]*kappa*cos(baseradians[radian]));

      radians = radians + (exp(-kc)*boost::math::cyl_bessel_i(0, kc) /
        (2*PI*exp(-kappa)*boost::math::cyl_bessel_i(0, kappa) *  
          exp(-pars[2])*boost::math::cyl_bessel_i(0, pars[2]))) *
          exp(kc - (kappa + pars[2]));
      
    }
    
    out[radian] = radians;
    radians = 0;
  }
  
  return(out/sum(out));
}
