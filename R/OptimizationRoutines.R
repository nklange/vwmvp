
one_two_pi <- 1/(2*pi)

ll_vp_sim <- function(pars, model, error_list, set_sizes, nsim, ...){
  
  precision <- pars[1]/(set_sizes^pars[2])
  parscont <- c(pars[3])
  
  if (grepl("RNplus",model)){
    parscont <- c(parscont,pars[4])
  }
  
  ll_fun <- match.fun(paste0("sim_fun_",model))
  
  out <- vector("numeric", length(set_sizes))
  
  for (i in seq_along(error_list)) {
    
    temp <- ll_fun(x = nsim, pars = c(precision[i],parscont),base_radians = base_radians)
    out[[i]] <- sum(log(temp[error_list[[i]]]))
  }
  
  return(sum(out))
}


ll_vp_numint <- function(pars, model, error_list, set_sizes){
  
  
  if (model %in% c("MK_F_RNplus","MK_F_RNminus")){
    
    #final value in K_range is non-integer K for precision(K,mKappa,alpha) for cases of K < SetSize
    #for K >= SetSize precision(SetSize,mKappa,alpha)
    SzK <- unique(c(set_sizes,pars[length(pars)]))
    K_range <- c(SzK[SzK <= pmin(max(SzK), pars[length(pars)])]) # if K > max(Sz) length of precision is higher than sz
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3:length(pars)])
    
  } else if (model %in% c("MK_FM_RNplus","MK_FM_RNminus","MK_FM2_RNplus","MK_FM2_RNminus")){
    
    #final value in K_range is integer ceiling(K) for precision(ceiling(K),mKappa,alpha) for cases of K < SetSize
    #second to final value in K_range is integer floor(K) for precision(floor(K),mKappa,alpha) for cases of K < SetSize
    # remaining values are set sizes smaller K: for K >= SetSize for precision(SetSize,mKappa,alpha)
    
    
    K_range <- c(set_sizes[set_sizes <= pmin(max(set_sizes), pars[length(pars)])], 
                 floor(pars[length(pars)]), ceiling(pars[length(pars)]))
    
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3:length(pars)])
    
    # mixture weights for floor(K),ceiling(K) from real part of non-integer K
    realK <- pars[length(pars)] - floor(pars[length(pars)])
    fmW <- c(1-realK,realK)
    
  } else if (model %in% c("MK_P_RNplus","MK_P_RNminus")){
    
    
    K_range <- c(0:max(set_sizes))
    poissW <- dpois(c(0:(max(K_range)-1)), pars[length(pars)], log = FALSE)
    poissW <- c(poissW,1-sum(poissW))
    
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3:(length(pars) - 1)])
    
  } else if (model %in% c("MK_U_RNplus","MK_U_RNminus")) {
    
    K_range <- c(0:max(set_sizes))
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3:length(pars)])
    
  } else if (model %in%  c("EP_RNplus","EP_RNminus")) {
    
    K_range <- set_sizes
    precision <- pars[1]/(K_range^pars[2])
    
    if (grepl("RNplus",model)) {
      parscont <- c(pars[3])
    } else {
      parscont <- c()
    }
    
  } else if (model %in% c("EP_F_RNplus","EP_F_RNminus")){
    
    #final value in K_range is non-integer K for precision(K,mKappa,alpha) for cases of K < SetSize
    #for K >= SetSize precision(SetSize,mKappa,alpha)
    SzK <- unique(c(set_sizes,pars[length(pars)]))
    K_range <- c(SzK[SzK <= pmin(max(SzK), pars[length(pars)])]) # if K > max(Sz) length of precision is higher than sz
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3:length(pars)])
    
  } else if (model %in% c("EP_FM_RNplus","EP_FM_RNminus")){
    
    K_range <- c(set_sizes[set_sizes <= pmin(max(set_sizes), pars[length(pars)])], 
                 floor(pars[length(pars)]), ceiling(pars[length(pars)]))
    
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3:length(pars)])
    
    # mixture weights for floor(K),ceiling(K) from real part of non-integer K
    realK <- pars[length(pars)] - floor(pars[length(pars)])
    fmW <- c(1-realK,realK)
    
  
  } else if (model %in% c("EP_P_RNplus","EP_P_RNminus")){
    
    
    K_range <- c(0:max(set_sizes))
    poissW <- dpois(c(0:(max(K_range)-1)), pars[length(pars)], log = FALSE)
    poissW <- c(poissW,1-sum(poissW))
    
    precision <- pars[1]/(K_range^pars[2])
    
    if (grepl("RNplus",model)) {
      parscont <- c(pars[3:(length(pars)-1)])
    } else {
      parscont <- c()
    }
    
    
  } else if (model %in% c("EP_U_RNplus","EP_U_RNminus")) {
    
    K_range <- c(0:max(set_sizes))
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3:length(pars)])
    
  } else if (model %in% c("SA_RNplus")){
    
    parscont <- pars[1]
    
  } else if (model %in% c("SA_F_RNplus","SA_U_RNplus","SA_F_RNminus","SA_U_RNminus")){
    
    precision <- pars[1]
    parscont <- c(pars[2:length(pars)])
    
  }  else if (model %in% c("SA_P_RNplus","SA_P_RNminus")){
    
    K_range <- c(0:max(set_sizes))
    poissW <- dpois(c(0:(max(K_range)-1)), pars[length(pars)], log = FALSE)
    poissW <- c(poissW,1-sum(poissW))
    
    precision <- pars[1]
    
    if (grepl("RNplus",model)){
      parscont <- c(pars[length(pars)-1])
    } else {
      parscont <- c()
    }
    
  } else if (model %in% c("MK_RNplus","MK_RNminus","J_RNplus","J_RNminus")) {
    
    K_range <- set_sizes
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3:length(pars)])
    
    
  }
  
  
  out <- vector("numeric", length(set_sizes))
  
  for (i in seq_along(error_list)) {
    
    if (model %in% c("MK_F_RNplus","MK_F_RNminus")){
      
      out[[i]] <- vp_f_routine(precision = precision, parscont = parscont, errors = error_list[[i]], 
                               set_sizes = set_sizes, sz = i, model = model)
      
    } else if (model %in% c("MK_FM_RNplus","MK_FM_RNminus")){
      
      out[[i]] <- vp_fm_routine(precision = precision, parscont = parscont, weights = fmW, errors = error_list[[i]],
                                set_sizes = set_sizes, sz = i, model = model)
      
    } else if (model %in% c("MK_FM2_RNplus","MK_FM2_RNminus")){
      
      out[[i]] <- vp_fm2_routine(precision = precision, parscont = parscont, weights = fmW, errors = error_list[[i]],
                                set_sizes = set_sizes, sz = i, model = model)
      
    } else if (model %in% c("MK_P_RNplus","MK_P_RNminus")) {
      
      out[[i]] <- vp_p_routine(precision = precision, parscont = parscont, weights = poissW, errors = error_list[[i]],
                               set_sizes = set_sizes, sz = i, model = model)
      
    } else if (model %in% c("MK_U_RNplus","MK_U_RNminus")) {
      
      out[[i]] <- vp_u_routine(precision = precision, parscont = parscont, errors = error_list[[i]],
                               set_sizes = set_sizes, sz = i,model = model)
      
    } else if (model %in% c("EP_RNplus","EP_RNminus")){
      
      out[[i]] <- ep_routine(pars = c(precision[i],parscont), errors = error_list[[i]],
                             model = model)
      
      
    } else if (model %in% c("EP_F_RNplus","EP_F_RNminus")){
      
      out[[i]] <- ep_f_routine(precision = precision, parscont = parscont, errors = error_list[[i]], 
                               set_sizes = set_sizes, sz = i, model = model)
      
      
    }  else if (model %in% c("EP_FM_RNplus","EP_FM_RNminus")){
      
      out[[i]] <- ep_fm_routine(precision = precision, parscont = parscont, weights = fmW, errors = error_list[[i]],
                                set_sizes = set_sizes, sz = i, model = model)
    
      } else if (model %in% c("EP_P_RNplus","EP_P_RNminus")) {
      
      out[[i]] <- ep_p_routine(precision = precision, parscont = parscont, weights = poissW, errors = error_list[[i]],
                               set_sizes = set_sizes, sz = i, model = model)
      
    } else if (model %in% c("EP_U_RNplus","EP_U_RNminus")) {
      
      out[[i]] <- ep_u_routine(precision = precision, parscont = parscont, errors = error_list[[i]],
                               set_sizes = set_sizes, sz = i, model = model)
      
    } else if (model %in% c("SA_RNplus")){
      
      out[[i]] <- sa_routine(pars = c(parscont), errors = error_list[[i]])
      
      
    } else if (model %in% c("SA_F_RNplus","SA_F_RNminus")){
      
      out[[i]] <- sa_f_routine(pars = c(precision,parscont), errors = error_list[[i]], 
                               set_sizes = set_sizes, sz = i, model = model)
      
      
    } else if (model %in% c("SA_U_RNplus","SA_U_RNminus")){
      
      out[[i]] <- sa_u_routine(pars = c(precision,parscont), errors = error_list[[i]], 
                               set_sizes = set_sizes, sz = i, model = model)
      
      
    } else if (model %in% c("SA_P_RNplus","SA_P_RNminus")){
      
      out[[i]] <- sa_p_routine(pars = c(precision,parscont), errors = error_list[[i]], 
                               weights = poissW, set_sizes = set_sizes, sz = i, model = model)
      
      
    } else  if (model %in% c("MK_RNplus","MK_RNminus","J_RNplus","J_RNminus")){
      
      out[[i]] <- vp_routine(precision = precision[i], parscont = parscont, errors = error_list[[i]],
                             model = model)
      
    }
    
    
  }
  
  return(sum(out))
}

vp_integration <- function(error, pars, coreFunction) {
  tryCatch(
    integrate(match.fun(coreFunction), 
              0, Inf, 
              pars,
              radian = error, stop.on.error = FALSE)$value, 
    error = function(e) tryCatch(
      integrate(match.fun(coreFunction), 
                0, Inf, 
                pars,
                radian = if (error == 0) {
                  circular(.Machine$double.xmin)  
                } else error, stop.on.error = FALSE)$value), 
    error = function(e) tryCatch(
      integrate(match.fun(coreFunction), 
                0, Inf,   
                pars,
                radian = if (error == 0) {
                  circular(.Machine$double.eps^2) 
                } else error, stop.on.error = FALSE)$value), error = function(e) NA)
}
rnplus_integration <- function(pars, errors){
  
  out <- vector("numeric",length(errors))
  kc <- vector("numeric",length(errors))
  
  kc <- sqrt(pars[2]^2 + pars[1]^2 + 2 * pars[2]*pars[1]*cos(errors))
  
  out <- 
    (besselI(kc,0,expon.scaled = TRUE) / 
       (2*pi*besselI(pars[1],0,expon.scaled = TRUE) * 
          besselI(pars[2],0,expon.scaled = TRUE))) *
    exp(kc - (pars[1] + pars[2]))
  
  return(out)
  
}
rnminus_integration <- function(pars, errors){
  
  out <- vector("numeric",length(errors))
  
  out <-   1 / (2 * pi * besselI(x = pars[1],nu = 0,expon.scaled = TRUE)) *
    (exp(cos(errors) - 1)) ^ pars[1]
  
  return(out)
}

uniformweights <- function(set_sizes, sz, parsK){
  
  ## calculate proportion of real part of K from total K
  unifw_real <- (parsK - floor(parsK)) / parsK
  ## calculate weight for integer K values from 0 to floor(K)
  unifw_int <- rep((1 - unifw_real)/ceiling(parsK), ceiling(parsK))
  ## weights for 0:ceiling(K)
  unifw_all <- c(unifw_int, unifw_real)
  ## restrict weights to current set_size and give remainder to last set_size
  unifw <- unifw_all[1:min(ceiling(parsK)+1,set_sizes[[sz]]+1)]
  unifw[length(unifw)] <- 1-sum(unifw[1:(length(unifw)-1)])
  return(unifw)
}

ep_routine <- function(pars,errors, model){
  
  
  modeltype <- regmatches(model, regexpr("[^_]+$", model))
  coreFunction <- paste0(tolower(modeltype),"_integration")
  
  
  out <- match.fun(coreFunction)(pars,errors)
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
}
ep_f_routine <- function(precision, parscont, errors, set_sizes, sz, model){
  
  
  modeltype <- regmatches(model, regexpr("[^_]+$", model))
  coreFunction <- paste0(tolower(modeltype),"_integration")
  
  
  pEncode <-  min(parscont[length(parscont)],set_sizes[[sz]]) / set_sizes[[sz]]
  
  if (modeltype == "RNplus") {
    pars <- c(precision[[pmin(sz, length(precision))]],parscont[c(1)])
  } else {
    pars <- c(precision[[pmin(sz, length(precision))]])
    
  }
  
  err <-match.fun(coreFunction)(pars,errors)
  
  out <- pEncode*err + (1-pEncode)*one_two_pi
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
}
ep_p_routine <- function(precision, parscont, weights, errors, set_sizes, sz, model){
  
  
  modeltype <- regmatches(model, regexpr("[^_]+$", model))
  coreFunction <- paste0(tolower(modeltype),"_integration")
  
  out <- vector("numeric", length(errors))
  
  out <- out + weights[[1]] * one_two_pi # for K = 0
  
  for (j in c(2:length(weights))){ # K = 1 to K = max(setsize) integer items
    
    K <- j - 1
    
    
    if (K < set_sizes[[sz]]){
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[j]],parscont)
      
      err <- match.fun(coreFunction)(pars,errors)
      
      pEncode <- K/set_sizes[[sz]]
      out <- out +  weights[[j]]*(pEncode*err + (1-pEncode)*one_two_pi)
      
    } else {
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[(set_sizes[[sz]] + 1)]],parscont)
      
      err <- match.fun(coreFunction)(pars,errors)
      
      out <- out +  weights[[j]]*err
    }
    
    
  }
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
}
ep_u_routine <- function(precision, parscont, errors, set_sizes, sz, model) {
  
  
  modeltype <- regmatches(model, regexpr("[^_]+$", model))
  coreFunction <- paste0(tolower(modeltype),"_integration")
  
  
  unifw <- uniformweights(set_sizes = set_sizes, sz = sz, parsK = parscont[length(parscont)])
  
  
  out <- vector("numeric", length(errors))
  
  
  out <- out + unifw[[1]] * one_two_pi # for K = 0
  
  for (j in c(2:length(unifw))){ # K = 1 to K = min(floor(Kpar),setsizes[[sz]]) integer items
    
    K <- j - 1
    
    
    if (K < set_sizes[[sz]]){
      
      err<- vector("numeric",length(errors))
      
      if(grepl("RNplus",model)){
        pars <- c(precision[[j]],parscont[c(1)])
      } else {
        pars <- precision[[j]]
      }
      
      err <- match.fun(coreFunction)(pars,errors)
      
      pEncode <- K/set_sizes[[sz]]
      out <- out +  unifw[[j]]*(pEncode*err + (1-pEncode)*one_two_pi)
      
    } else {
      
      err<- vector("numeric",length(errors))
      
      if(grepl("RNplus",model)){
        pars <- c(precision[[(set_sizes[[sz]] +1)]],parscont[c(1)])
      } else {
        pars <- precision[[(set_sizes[[sz]] +1)]]
      }
      
      err <- match.fun(coreFunction)(pars,errors)
      
      out <- out +  unifw[[j]]*err
    }
    
    
  }
  
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
  
}
ep_fm_routine <- function(precision, parscont, weights, errors, set_sizes, sz, model){
  
  modeltype <- regmatches(model, regexpr("[^_]+$", model))
  coreFunction <- paste0(tolower(modeltype),"_integration")
  
  out <- vector("numeric", length(errors))
  Ss <- set_sizes[[sz]]
  pEncode <-  min(parscont[length(parscont)],set_sizes[[sz]]) / set_sizes[[sz]]
  condpEncode<- c(floor(parscont[length(parscont)])/Ss,ceiling(parscont[length(parscont)])/Ss)
  
  
  if (pEncode < 1){
    
    for (FMind in seq_along(weights)){
      
      err<- vector("numeric",length(errors))
      
      pars <- c(tail(precision,2)[[FMind]],parscont[c(1:(length(parscont) - 1))])
      
      err <- match.fun(coreFunction)(pars,errors)
      
      out <-  out + weights[[FMind]] * (condpEncode[[FMind]]*err + (1-condpEncode[[FMind]])*one_two_pi)
    }
    
    
  } else {
    
    
    err<- vector("numeric",length(errors))
    pars <- c(precision[[sz]],parscont[c(1:(length(parscont) - 1))])
    
    err <- match.fun(coreFunction)(pars,errors)
    
    out <- err
    
  }
  
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
}

vp_routine <- function(precision, parscont, errors, model) {
  
  coreFunction <- paste0("cint_fun_",model)
  
  out <- vector("numeric", length(errors))
  
  pars <- c(precision,parscont)
  
  
  for (i in seq_along(out)) {
    out[i] <- vp_integration(error = errors[i], pars = pars, 
                             coreFunction = coreFunction)
  }
  
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
  
}
vp_u_routine <- function(precision, parscont, errors, set_sizes, sz, model) {
  
  modeltype <- regmatches(model, regexpr("_[^_]+$", model))
  coreFunction <- paste0("cint_fun_MK",modeltype)
  
  unifw <- uniformweights(set_sizes = set_sizes, sz = sz, parsK = parscont[length(parscont)])
  
  out <- vector("numeric", length(errors))
  for (j in c(1:length(unifw))){ # K = 0 to K = min(floor(Kpar),setsizes[[sz]]) integer items
    
    K <- j - 1
    
    
    if (K < set_sizes[[sz]]){
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[j]],parscont[c(1:(length(parscont)-1))])
      
      for (i in seq_along(err)) {
        
        err[i] <- vp_integration(error = errors[i], pars = pars, 
                                 coreFunction = coreFunction)
      }
      
      pEncode <- K/set_sizes[[sz]]
      out <- out +  unifw[[j]]*(pEncode*err + (1-pEncode)*one_two_pi)
      
    } else {
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[(set_sizes[[sz]] +1)]],parscont[c(1:(length(parscont)-1))])
      
      for (i in seq_along(err)) {
        
        err[i] <- vp_integration(error = errors[i], pars = pars, 
                                 coreFunction = coreFunction)
      }
      
      out <- out +  unifw[[j]]*err
    }
    
    
  }
  
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
  
}
vp_p_routine <- function(precision, parscont, weights, errors, set_sizes, sz, model) {
  
  
  modeltype <- regmatches(model, regexpr("_[^_]+$", model))
  coreFunction <- paste0("cint_fun_MK",modeltype)
  
  out <- vector("numeric", length(errors))
  for (j in c(1:length(weights))){ # K = 0 to K = max(setsize) integer items
    
    K <- j - 1 #to start from Kitems = 0
    
    
    if (K < set_sizes[[sz]]){
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[j]],parscont)
      
      for (i in seq_along(err)) {
        
        err[i] <- vp_integration(error = errors[i], pars = pars, 
                                 coreFunction = coreFunction)
      }
      
      pEncode <- K/set_sizes[[sz]]
      out <- out +  weights[[j]]*(pEncode*err + (1-pEncode)*one_two_pi)
      
    } else {
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[(set_sizes[[sz]] + 1)]],parscont)
      
      for (i in seq_along(err)) {
        
        err[i] <- vp_integration(error = errors[i], pars = pars, 
                                 coreFunction = coreFunction)
      }
      
      out <- out +  weights[[j]]*err
    }
    
    
  }
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
  
}
vp_f_routine <- function(precision, parscont, errors, set_sizes, sz, model) {
  
  modeltype <- regmatches(model, regexpr("_[^_]+$", model))
  coreFunction <- paste0("cint_fun_MK",modeltype)
  
  pEncode <-  min(parscont[length(parscont)],set_sizes[[sz]]) / set_sizes[[sz]]
  pars <- c(precision[[pmin(sz, length(precision))]],parscont[c(seq(1,length(parscont)-1,1))])
  
  out <- vector("numeric", length(errors))
  err<- vector("numeric",length(errors))
  
  
  for (i in seq_along(err)) {
    err[i] <- vp_integration(error = errors[i], pars = pars, 
                             coreFunction = coreFunction)
  }
  
  
  out <- pEncode*err + (1-pEncode)*one_two_pi
  
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
  
}
vp_fm_routine <- function(precision, parscont, weights, errors, set_sizes, sz, model) {
  
  modeltype <- regmatches(model, regexpr("_[^_]+$", model))
  coreFunction <- paste0("cint_fun_MK",modeltype)
  
  out <- vector("numeric", length(errors))
  Ss <- set_sizes[[sz]]
  pEncode <-  min(parscont[length(parscont)],Ss) / Ss
  condpEncode<- c(floor(parscont[length(parscont)])/Ss,ceiling(parscont[length(parscont)])/Ss)
  
  if (pEncode < 1){
    
    for (FMind in seq_along(weights)){
      
      err<- vector("numeric",length(errors))
      pars <- c(tail(precision,2)[[FMind]],parscont[c(1:(length(parscont) - 1))])
      
      for (i in seq_along(err)) {
        err[i] <- vp_integration(error = errors[i], pars = pars, 
                                 coreFunction = coreFunction)
        
      }
      out <- out + weights[[FMind]] * (condpEncode[[FMind]]*err + (1-condpEncode[[FMind]])*one_two_pi)
    }
    
    
  } else {
    
    
    err<- vector("numeric",length(errors))
    pars <- c(precision[[sz]],parscont[c(1:(length(parscont) - 1))])
    
    for (i in seq_along(err)) {
      err[i] <- vp_integration(error = errors[i], pars = pars, 
                               coreFunction = coreFunction)
      
    }
    out <- err
    
  }
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
  
}
vp_fm2_routine <- function(precision, parscont, weights, errors, set_sizes, sz, model) {
  
  modeltype <- regmatches(model, regexpr("_[^_]+$", model))
  coreFunction <- paste0("cint_fun_MK",modeltype)
  
  out <- vector("numeric", length(errors))
  outweighed <- out
  pEncode <-  min(parscont[length(parscont)],set_sizes[[sz]]) / set_sizes[[sz]]
  
  if (pEncode < 1){
    
    for (FMind in seq_along(weights)){
      
      err<- vector("numeric",length(errors))
      pars <- c(tail(precision,2)[[FMind]],parscont[c(1:(length(parscont) - 1))])
      
      for (i in seq_along(err)) {
        err[i] <- vp_integration(error = errors[i], pars = pars, 
                                 coreFunction = coreFunction)
        
      }
      outweighed <- outweighed + weights[[FMind]] * err
    }
    
    out <- pEncode * outweighed + (1-pEncode)*one_two_pi
    
  } else {
    
    
    err<- vector("numeric",length(errors))
    pars <- c(precision[[sz]],parscont[c(1:(length(parscont) - 1))])
    
    for (i in seq_along(err)) {
      err[i] <- vp_integration(error = errors[i], pars = pars, 
                               coreFunction = coreFunction)
      
    }
    out <- err
    
  }
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
}

sa_routine <- function(pars, errors){
  
  out <- rnminus_integration(pars,errors)
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
}
sa_f_routine <- function(pars, errors, set_sizes, sz, model){
  
  modeltype <- regmatches(model, regexpr("[^_]+$", model))
  coreFunction <- paste0(tolower(modeltype),"_integration")
  
  Sz <- set_sizes[[sz]]
  pEncode <-  pars[length(pars)] / Sz
  
  if (pEncode <= 1){
    
    pars_lim <- pars[c(1:(length(pars) - 1))]
    
    err <- match.fun(coreFunction)(pars = pars_lim,errors)
    
    out <- pEncode*err + (1-pEncode)*one_two_pi
    
  } else {
    
    phigh <- (pars[length(pars)] %% Sz)/Sz # proportion of items that are encoded with the higher precision
    kappa_low <- pars[1] * floor(pEncode) # precision for 1 item times the increase in precision given by K exceeding Sz
    kappa_high <- pars[1] * (floor(pEncode)+1) # as above but extra precision if single quanta of precision don't neatly map onto items
    
    if (modeltype == "RNplus"){
      parshigh <- c(kappa_high,pars[2])
      parslow <- c(kappa_low,pars[2])
    } else {
      parshigh <- kappa_high
      parslow <- kappa_low
    }
    err_high <- match.fun(coreFunction)(pars = parshigh, errors)
    err_low <- match.fun(coreFunction)(pars = parslow, errors)
    
    out <- phigh * err_high + (1 - phigh) * err_low
    
  }
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
}
sa_u_routine <- function(pars, errors, set_sizes, sz, model){
  
  modeltype <- regmatches(model, regexpr("[^_]+$", model))
  coreFunction <- paste0(tolower(modeltype),"_integration")
  
  
  unifw <- uniformweights(set_sizes = set_sizes, sz = sz, parsK = pars[length(pars)])
  
  
  out <- vector("numeric", length(errors))
  
  
  out <- out + unifw[[1]] * one_two_pi # for K = 0
  
  
  for (j in c(2:length(unifw))){ # K = 1 to K = min(floor(Kpar),setsizes[[sz]]) integer items
    
    K <- j - 1
    
    Sz <- set_sizes[[sz]]
    pEncode <-  K / Sz
    
    if (pEncode <= 1){
      
      pars_lim <- pars[c(1:(length(pars)-1))]
      
      err <- match.fun(coreFunction)(pars = pars_lim,errors)
      
      out <- out + unifw[[j]] * (pEncode*err + (1-pEncode)*one_two_pi)
      
    } else {
      
      phigh <- (K %% Sz)/Sz # proportion of items that are encoded with the higher precision
      kappa_low <- pars[1] * floor(pEncode) # precision for 1 item times the increase in precision given by K exceeding Sz
      kappa_high <- pars[1] * (floor(pEncode)+1) # as above but extra precision if single quanta of precision don't neatly map onto items
      
      if (grepl("RNplus",model)) {
        pars_high <- c(kappa_high,pars[2])
        pars_low <- c(kappa_low,pars[2])
      } else {
        pars_high <- kappa_high
        pars_low <- kappa_low
      }
      
      err_high <- match.fun(coreFunction)(pars = pars_high, errors)
      err_low <- match.fun(coreFunction)(pars = pars_low, errors)
      
      out <- out + unifw[[j]] * (phigh * err_high + (1 - phigh) * err_low)
    }
    
  }
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
}
sa_p_routine <- function(pars, errors, weights, set_sizes, sz, model){
  
  modeltype <- regmatches(model, regexpr("[^_]+$", model))
  coreFunction <- paste0(tolower(modeltype),"_integration")
  
  out <- vector("numeric", length(errors))
  
  
  out <- out + weights[[1]] * one_two_pi # for K = 0
  
  
  for (j in c(2:length(weights))){ # K = 1 to K = min(floor(Kpar),setsizes[[sz]]) integer items
    
    K <- j - 1
    
    Sz <- set_sizes[[sz]]
    pEncode <-  K / Sz
    
    if (pEncode <= 1){
      
      err <- match.fun(coreFunction)(pars = pars, errors)
      
      out <- out + weights[[j]] * (pEncode*err + (1-pEncode)*one_two_pi)
      
    } else {
      
      phigh <- (K %% Sz)/Sz # proportion of items that are encoded with the higher precision
      kappa_low <- pars[1] * floor(pEncode) # precision for 1 item times the increase in precision given by K exceeding Sz
      kappa_high <- pars[1] * (floor(pEncode)+1) # as above but extra precision if single quanta of precision don't neatly map onto items
      
      if (grepl("RNplus",model)) {
        pars_high <- c(kappa_high,pars[2])
        pars_low <- c(kappa_low,pars[2])
      } else {
        pars_high <- kappa_high
        pars_low <- kappa_low
      }
      
      err_high <- match.fun(coreFunction)(pars = pars_high, errors)
      err_low <- match.fun(coreFunction)(pars = pars_low, errors)
      
      out <- out + weights[[j]] * (phigh * err_high + (1 - phigh) * err_low)
    }
    
  }
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
}
