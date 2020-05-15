predict_data <- function(data, model, pars){
  
  dp <- prep_data(data)
  error_list <- dp$datalist
  set_sizes <- dp$set_sizes
  
  # Note: parameters are not named in function, so order has to be
  ## mkappa1/j1bar, alpha, tau, kappa_r, K
  ## parameters can be left off if appropriate, i.e.
  ### MK_F_RNplus: mkappa1/j1bar, alpha, tau, kappa_r, K
  ### EP_F_RNplus: mkappa1/j1bar, alpha, kappa_r, K
  ### SA_F_RNplus: mkappa1/j1bar, kappa_r, K
  
  param <- prep_parameters(pars, model, set_sizes)
  
  precision <- param[[1]]
  parscont <- param[[2]]
  weights <- param[[3]]
  
  out_list <- vector("list", length(set_sizes))
  
  for (i in seq_along(error_list)) {
    
    if (model %in% c("MK_F_RNplus","MK_F_RNminus")){
      
      out <- vp_f_routine(precision = precision, parscont = parscont, errors = error_list[[i]], 
                               set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
    } else if (model %in% c("MK_FM_RNplus","MK_FM_RNminus")){
      
      out <- vp_fm_routine(precision = precision, parscont = parscont, weights = weights, errors = error_list[[i]],
                                set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
    } else if (model %in% c("MK_FM2_RNplus","MK_FM2_RNminus")){
      
      out <- vp_fm2_routine(precision = precision, parscont = parscont, weights = weights, errors = error_list[[i]],
                                 set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
    } else if (model %in% c("MK_P_RNplus","MK_P_RNminus")) {
      
      out <- vp_p_routine(precision = precision, parscont = parscont[-length(parscont)], weights = weights, errors = error_list[[i]],
                               set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
    } else if (model %in% c("MK_P2_RNplus","MK_P2_RNminus")) {
      
      out <- vp_p2_routine(precision = precision, parscont = parscont, weights = weights, errors = error_list[[i]],
                                set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
    } else if (model %in% c("MK_U_RNplus","MK_U_RNminus")) {
      
      out <- vp_u_routine(precision = precision, parscont = parscont, errors = error_list[[i]],
                               set_sizes = set_sizes, sz = i,model = model, predictOrLL = "predict")
      
    } else if (model %in% c("MK_U2_RNplus","MK_U2_RNminus")) {
      
      out <- vp_u2_routine(precision = precision, parscont = parscont, errors = error_list[[i]],
                                set_sizes = set_sizes, sz = i,model = model, predictOrLL = "predict")
      
    } else if (model %in% c("EP_RNplus","EP_RNminus")){
      
      out <- ep_routine(pars = c(precision[i],parscont), errors = error_list[[i]],
                             model = model, predictOrLL = "predict")
      
      
    } else if (model %in% c("EP_F_RNplus","EP_F_RNminus")){
      
      out <- ep_f_routine(precision = precision, parscont = parscont, errors = error_list[[i]], 
                               set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
      
    }  else if (model %in% c("EP_FM_RNplus","EP_FM_RNminus")){
      
      out <- ep_fm_routine(precision = precision, parscont = parscont, weights = weights, errors = error_list[[i]],
                                set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
    }  else if (model %in% c("EP_FM2_RNplus","EP_FM2_RNminus")){
      
      out<- ep_fm2_routine(precision = precision, parscont = parscont, weights = weights, errors = error_list[[i]],
                                 set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
    } else if (model %in% c("EP_P_RNplus","EP_P_RNminus")) {
      
      # don't need realK in EP_P models, only in EP_P2 models
      if (model == "EP_P_RNplus"){
        parscont <- parscont[1:(length(parscont)-1)]
      } else if (model == "EP_P_RNminus"){
        parscont <- c()
      }
      
      out <- ep_p_routine(precision = precision, parscont = parscont, weights = weights, errors = error_list[[i]],
                               set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
    } else if (model %in% c("EP_P2_RNplus","EP_P2_RNminus")) {
      
      out <- ep_p2_routine(precision = precision, parscont = parscont, weights = weights, errors = error_list[[i]],
                                set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
    }  else if (model %in% c("EP_U_RNplus","EP_U_RNminus")) {
      
      out <- ep_u_routine(precision = precision, parscont = parscont, errors = error_list[[i]],
                               set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
    } else if (model %in% c("EP_U2_RNplus","EP_U2_RNminus")) {
      
      out <- ep_u2_routine(precision = precision, parscont = parscont, errors = error_list[[i]],
                                set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
    } else if (model %in% c("SA_RNplus")){
      
      out <- sa_routine(pars = c(parscont), errors = error_list[[i]], predictOrLL = "predict")
      
      
    } else if (model %in% c("SA_F_RNplus","SA_F_RNminus")){
      
      out <- sa_f_routine(pars = c(precision,parscont), errors = error_list[[i]], 
                               set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
      
    } else if (model %in% c("SA_U_RNplus","SA_U_RNminus")){
      
      out <- sa_u_routine(pars = c(precision,parscont), errors = error_list[[i]], 
                               set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
      
    } else if (model %in% c("SA_P_RNplus","SA_P_RNminus")){
      
      out <- sa_p_routine(pars = c(precision,parscont), errors = error_list[[i]], 
                               weights = weights, set_sizes = set_sizes, sz = i, model = model, predictOrLL = "predict")
      
      
    } else  if (model %in% c("MK_RNplus","MK_RNminus","J_RNplus","J_RNminus")){
      
      out <- vp_routine(precision = precision[i], parscont = parscont, errors = error_list[[i]],
                             model = model, predictOrLL = "predict")
      
    } else if (model %in% c("UVM")){
      
      out <- uvm_routine(pars = parscont, errors = error_list[[i]], model = model, predictOrLL = "predict")
    }
    
    
    out_list[[i]] <- bind_cols(tibble::tibble(prediction = out),
                tibble::tibble(data = as.numeric(error_list[[i]])),
                tibble::tibble(setsize = rep(set_sizes[[i]],length(out))),
                tibble::tibble(model = rep(model,length(out))),
                tibble::tibble(id = rep(dp$id,length(out))))
      
      
  }
  
  
  dplyr::bind_rows(out_list)
}

generate_data <- function(model, par, trials, set_sizes){
  
  if (grepl("MK",model) | grepl("J",model)){
    Err <- generate_VP_data(trials = trials, set_sizes = set_sizes, model = model, par = par)
  } else if (grepl("EP",model)){
    Err <- generate_EP_data(trials = trials, set_sizes = set_sizes, model = model, par = par)
  } else if (grepl("SA",model)){
    Err <- generate_SA_data(trials = trials, set_sizes = set_sizes, model = model, par = par)
  }
  return(Err)
}

generate_VP_data <- function(trials, set_sizes, model, par){
  
  Err <- NULL
  
  # Name Parameters
  meanprec <- par[[1]]
  alpha <- par[[2]]
  tau <- par[[3]]
  if (grepl("RNplus",model)) {
    kappa_r <- par[[4]]
  }
  if (model %in% c("MK_FM_RNplus","MK_U_RNplus","MK_P_RNplus")){
    Kpar <- par[[5]]
  } 
  
  if (model %in% c("MK_FM_RNminus","MK_U_RNminus","MK_P_RNminus")){
    Kpar <- par[[4]]
  } 
  
  
  for (setsizeInd in seq_along(set_sizes)){
    trialsSs <- trials[[setsizeInd]]
    
    # draw target stimuli for all trials
    StimuliTrial <- runif(trialsSs)*2*pi-pi
    
    
    # Get type of guessing right
    if (model %in% c("MK_FM_RNminus","MK_FM_RNplus")){
      KTrial <- ifelse(runif(trialsSs) < Kpar- floor(Kpar),ceiling(Kpar),floor(Kpar))
      
    } else if (model %in% c("MK_P_RNminus","MK_P_RNplus")){
      KTrial <- rpois(trialsSs,Kpar)
      
    } else if (model %in% c("MK_U_RNminus","MK_U_RNplus")){
      
      KTrial <- round(runif(trialsSs,min = 0, max = Kpar))
      
    } else {
      
      KTrial <- rep(Inf,trialsSs)
      
    }
    
    
    Ss <- set_sizes[[setsizeInd]]
    
    # determine precision for trials
    
    Errors <- vector("numeric",length(StimuliTrial))
    
    for (i in seq_along(StimuliTrial)){
      
      kappabar <-  meanprec/(min(KTrial[i],Ss)^(alpha))

      
      memory <- ifelse(KTrial[i] > Ss,TRUE,ifelse(runif(1) < KTrial[i]/Ss,TRUE,FALSE))
      # response from memory if K > Setsize or, with probability K/N, i.e. if runif is smaller than K/N
      
      if (memory){
        
        # VP models
        # generate kappa from gamma for each trial
        kappagam <- rgamma(1,shape = kappabar/tau,scale = tau)
        
        if (model %in% c("J_RNminus","J_RNplus")){
          kappa <- KappafromJ(kappagam)
        } else {
          kappa <- kappagam
        }
        
        
        if (kappa < 1e-6) {
          
          # make a guess because rvonmises breaks, and it's basically guessing
          sbar <- circular::circular(runif(1)*2*pi-pi)
          
        } else {
          
          # remember stimulus with some precision kappa
          sbar <- circular::rvonmises(1,mu = circular::circular(StimuliTrial[[i]]),kappa,control.circular=list(units="radians"))
         
        }
      } else {
        sbar <- circular::circular(runif(1)*2*pi-pi)
      }
      
      if(grepl("RNplus",model)){
        # add response noise kappa_r to memory
        # calculate error as distance from presented target
        err <- StimuliTrial[[i]] - circular::rvonmises(1,mu=sbar,kappa_r)
      } else {
        err <- StimuliTrial[[i]] - sbar
      }
      
      # make sure distance is represented on a circle
      Errors[[i]] <- ifelse(err < -pi, err + 2*pi, 
                            ifelse(err > pi, err- 2*pi, err))
      
    } 
    
    SzErrors <- bind_cols(tibble::tibble(error_0 = Errors),
                          tibble::tibble(set_size = rep(set_sizes[[setsizeInd]],trialsSs))
    )
    Err <- bind_rows(Err,SzErrors)
  }
  
  
  return(Err)
}
