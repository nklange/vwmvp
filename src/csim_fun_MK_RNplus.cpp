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
  double kc;
  int lr = baseradians.size();
  NumericVector out(lr);
  double radians = 0;
  
  kappa = Rcpp::rgamma(x,pars[0]/pars[1], pars[1]);
  
  for (radian = 0; radian < lr; radian++){
    
    for (i = 0; i < x; i++) {

      kc = sqrt(pow(pars[2], 2) + pow(kappa[i], 2) + 2 * pars[2]*kappa[i]*cos(baseradians[radian]));
      
      // Rcout << kc << std::endl;
      
      radians += ((R::bessel_i(kc, 0, 2) /
        (2*PI*R::bessel_i(kappa[i], 0, 2) *  
          R::bessel_i(pars[2], 0, 2))) *
          exp(kc - (kappa[i] + pars[2])));
      
    }
    
    out[radian] = radians;
    radians = 0;
  }
  
  return(out/sum(out));
}
