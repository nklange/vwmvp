#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector csim_fun_MK_RNplus(
    double x,
    NumericVector pars, 
    NumericVector baseradians) {
  
  int i;
  int radian;
  NumericVector kappa;
  NumericVector kc;
  int lr = baseradians.size();
  NumericVector out(lr);
  double radians = 0;
  
  kappa = Rcpp::rgamma(x,pars[0]/pars[1], pars[1]);
  
  for (radian = 0; radian < lr; radian++){
    
    kc = sqrt(pow(pars[2], 2) + pow(kappa, 2) + 2 * pars[2]*kappa*cos(baseradians[radian]));
    
    for (i = 0; i < x; i++) {+
      
      // Rcout << kc << std::endl;
      
      radians += (R::bessel_i(kc[i], 0, 2) /
        (2*PI*R::bessel_i(kappa[i], 0, 2) *  
          R::bessel_i(pars[2], 0, 2))) *
          exp(kc[i] - (kappa[i] + pars[2]));
      
    }
    
    out[radian] = radians;
    radians = 0;
  }
  
  return((out/sum(out)) * 0.5);
}