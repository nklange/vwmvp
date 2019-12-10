// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
# include <cmath>
# include <complex>
# include <cstdlib>
# include <cstring>
# include <iomanip>
# include <iostream>

using namespace Rcpp;
using namespace std;

# include "fn.hpp"


// [[Rcpp::export]]
float r4_besi0e ( float x );

// [[Rcpp::export]]
NumericVector csim_fun_MK_RNplus(
    double x,
    NumericVector pars, 
    NumericVector baseradians) {
  
  int i;
  int radian;
  NumericVector kappa;
  double kc;
  int lr = baseradians.size();
  NumericVector out(lr);
  double radians = 0;
  
  kappa = Rcpp::rgamma(x,pars[0]/pars[1], pars[1]);
  
  for (radian = 0; radian < lr; radian++){
    
    for (i = 0; i < x; i++) {

      kc = sqrt(pow(pars[2], 2) + pow(kappa[i], 2) + 2 * pars[2]*kappa[i]*cos(baseradians[radian]));
      
      // Rcout << kc << std::endl;
      
      radians = radians + (exp(-kc)*boost::math::cyl_bessel_i(0, kc) /
        (2*PI*exp(-kappa[i])*boost::math::cyl_bessel_i(0, kappa[i]) *  
          exp(-pars[2])*boost::math::cyl_bessel_i(0, pars[2]))) *
          exp(kc - (kappa[i] + pars[2]));
      
    }
    
    out[radian] = radians;
    radians = 0;
  }
  
  return(out/sum(out));
}
