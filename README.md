# vwmvp

Collection of functions to fit (mostly) variable precision models to data produced by the delayed estimation paradigm in visual working memory. Wrapped as a package to make its use across various rstudio projects/files easier.

# Pre-requisites for data

Data needs to be a dataframe/tibble/or similar with at a minimum three columns, where each row represents one trial (long-format): column 'id': participant number, column 'set_size': set size of the trial, column 'error_0': deviance of participant response from target feature in radians, ranging from -pi to pi. Additional columns that are reproduced in the output are: exp (experiment id), cvid (separate id for crossvalidation folds), leftout (indicating fold for crossvalidation).

```
kappa <- c(20,10,2)
samplesize <- 100

SampleData <- function(kappa,samplesize){

  sample <- as.numeric(circular::rvonmises(samplesize,mu=circular::circular(0),kappa=kappa,control.circular = list(units = "radians")))
  sample_circls <- ifelse(sample < -pi,sample + 2*pi,ifelse(sample > pi, sample - 2*pi,sample))
  return(sample_circls)
}

# data format required for fitting, at a minimum:

data <- data.frame(id = "Test",
                   set_size = rep(c(1,2,6),each=samplesize),
                   error_0 = c(SampleData(kappa[[1]],samplesize),
                              SampleData(kappa[[2]],samplesize),
                              SampleData(kappa[[3]],samplesize))))
```

# Included models

Proper Term | Name in package
-----|--------------
VP(J)A-    | J_RNminus      
VP(J)A+        | J_RNplus
VP(&kappa;)A- | MK_RNminus
VP(&kappa;)A+ | MK_RNplus
VP(&kappa;)F- | MK_FM_RNminus
VP(&kappa;)F+ | MK_FM_RNplus
VP(&kappa;)P- | MK_P_RNminus
VP(&kappa;)P+ | MK_P_RNplus
VP(&kappa;)U- | MK_U_RNminus
VP(&kappa;)U+ | MK_U_RNplus
EP(&kappa;)A- | EP_RNminus
EP(&kappa;)A+ | EP_RNplus
EP(&kappa;)F- | EP_FM_RNminus
EP(&kappa;)F+ | EP_FM_RNplus
EP(&kappa;)P- | EP_P_RNminus
EP(&kappa;)P+ | EP_P_RNplus
EP(&kappa;)U- | EP_U_RNminus
EP(&kappa;)U+ | EP_U_RNplus
SA(&kappa;)A+ | SA_RNplus
SA(&kappa;)F- | SA_F_RNminus
SA(&kappa;)F+ | SA_F_RNplus
SA(&kappa;)P- | SA_P_RNminus
SA(&kappa;)P+ | SA_P_RNplus
SA(&kappa;)U- | SA_U_RNminus
SA(&kappa;)U+ | SA_U_RNplus

where:

* VP: variable precision
* EP: equal precision
* SA: slots and averaging
* &kappa;: mean memory precision parameterized directly as $\kappa$
* J: mean memory precision parameterized as Fisher information
* A: unlimited memory capacity
* F/FM: fixed memory capacity limit
* P: variable, poisson-distributed capacity limit
* U: variable, uniform-distributed capacity limit
* -: no response noise
* +: response noise

Additionally, some derivatives with slightly different implementation of the capacity limit for the VP models (F, F2, P2, U2) and a base VP model without set size effect, e.g., without the $\alpha$ parameter (VPnosetsize, VPplusnosetsize).

# Fitting routine

The top-level fitting function is **vwmvp::FitVP()** (in FittingAlgorithms.R). In principle, this can be used for the simulation approach for VP(&phi;)A&plusmn; models for the numerical integration approach for all models (where necessary). 

For the simulation approach (method = "sim"):

```
rep <- 20
seqrun <- 5
nsim <- 1500

vwmvp::FitVP(data = data, model = "MK_RNminus", rep=rep, method="sim", seqrun=seqrun, nsim=nsim)

```
where **rep** gives the number of model fitting runs, **seqrun** (default = 5) indicates the number of consecutive evaluations without improvement that will lead to the termination of the fitting run, and **nsim** (default = 1500) indicates the number of samples taken from the Gamma distribution.


For the numerical integration approach (method = "numint"):

```
rep <- 20

vwmvp::FitVP(data = data, model = "MK_RNminus", rep=rep, method = "numint")

```

The top-level function then internally draws on a number of different functions:

* **vwmvp::prep_data()** (in HelperFunctions.R): prepares data
* **vwmvp::get_start_vp()** (in HelperFunctions.R): generate starting parameters
* **vwmvp::fit_one_vp_nlminb()**/**vwmvwp::fit_one_vp_ga()** (in FittingAlgorithms.R): calls nlminb for numerical integration and GA::ga for simulation approach
* **vwmvp::ll_vp_numint** / **vwmvp::ll_vp_sim** (in OptimizationRoutines.R): hands down the model fitting to the model-specific routine for numerical integration and simulation approach respectively
* for VP models: **vwmvp::vp_routine** (and similar) that calls on **vwmvp::vp_integration** as the integration routine which in turn calls on the RCPP files for a standard VP model with/without response noise, for example, **vwmvp::cint_fun_MK_RNminus**
* additional functions for weights in models with limited capacity, and functions for non-VP models


Given data such as the one provided in the example above, the numerical integration approach will produce a tibble as below:

```

# A tibble: 3 x 17
   alpha mkappa1     tau objective convergence iterations message time  model id    leftout   rep exp   cvid  s_alpha s_mkappa1
   <dbl>   <dbl>   <dbl>     <dbl>       <int>      <int> <chr>   <drt> <chr> <fct>   <dbl> <int> <lgl> <fct>   <dbl>     <dbl>
1 1.2031  20.870 0.86974    164.65           0         14 relati~  9.4~ MK_R~ Test        0     1 NA    Test  1.4338     41.240
2 1.1994  20.863 0.88770    164.65           0         22 relati~ 16.0~ MK_R~ Test        0     2 NA    Test  1.1968     46.565
3 1.2030  20.871 0.87035    164.65           0         20 relati~ 12.8~ MK_R~ Test        0     3 NA    Test  0.81187    55.816
# ... with 1 more variable: s_tau <dbl>

```

where: 

* alpha, mkappa1,J,tau,kappa_r,K: parameter estimates associated with the fit
* s_alpha, s_mkappa1, s_J, s_tau, s_kappa_r, s_K: starting parameters
* objective: negative log likelihood
* convergence: convergence message from nlminb
* interations: number of evaluations of the likelihood in that run
* message: message associated with that fitting run
* time: elapsed time for that fitting run
* model: model fitted
* id: input id
* leftout: default = 0, unless specified 
* rep: number of run (max(rep) is **rep** specified in the top-level function)
* exp: default = NA, unless specified
* cvid: default = id, unless specified

## Interpreting K

The K parameter is estimated in models with limited memory capacity. The value produced as the parameter estimate can be between 0 and &infin;. In P/P2 models, K indicates the mean of the poisson distribution, formally Pois(&lambda; = K). In F/FM/FM2 models, it indicates the fixed capacity limit. In these models the limit cannot exceed the maximum set size in an experiment, i.e., as K >= max(set size), all items are remembered. For parameter estimates of K > max(set size) therefore K = &infin;. In the U/U2 models we discuss K as representing the mean of the uniform distribution (analogously to the poisson-distributed capacity limit), where the limit can range from 0 to 2K, with mean = K. Throughout the functions collected here (with one exception), we have actually implemented a range of 0 to K, with mean = K/2. When aiming to report the K parameter estimate, the given value of the parameter estimate should therefore be divided by 2 to represent the mean of the distribution.

# Making predictions

There are two approaches to make predictions using the functions.

In **vwmvp::generate_data()** (in PredictionRoutines.R), trials for a given model and parameters are simulated directly.

In **vwmvp::predict_data()** (in PredictionRoutines.R), the probability distribution associated with a model/parameters is reproduced using the fitting routines. By sampling from this, trials can be simulated.

[vwmvp_predictvsgenerate.R](https://github.com/nklange/InferenceModelComparison/blob/master/vwmvp_predictvsgenerate.R) shows that both approaches lead to approximately equivalent simulations, given a sufficiently high number of trials. Note: In vwmvp::generate_data(), the limits in U models is defined as ranging from 0 to 2K (versus 0 to K in vwmvp::predict_data()). Hence, for equivalent results, the input-K value in vwmvp::predict_data() has to be double that of the value in vwmvp::generate_data().