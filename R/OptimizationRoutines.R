
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
  

  
  if (model %in% c("MK_F_RNplus")){
    
    #K_range <- set_sizes[set_sizes <= pmin(max(set_sizes), pars[5])]
    
    #final value in K_range is non-integer K for precision(K,mKappa,alpha) for cases of K < SetSize
    #for K >= SetSize precision(SetSize,mKappa,alpha)
    K_range <- c(set_sizes[set_sizes <= pmin(max(set_sizes), pars[5])],pars[5])
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3],pars[4],pars[5])
    
  } else if (model %in% c("MK_FM_RNplus")){

    #final value in K_range is integer ceiling(K) for precision(ceiling(K),mKappa,alpha) for cases of K < SetSize
    #second to final value in K_range is integer floor(K) for precision(floor(K),mKappa,alpha) for cases of K < SetSize
    # remaining values are set sizes smaller K: for K >= SetSize for precision(SetSize,mKappa,alpha)
    K_range <- c(set_sizes[set_sizes <= pmin(max(set_sizes), pars[5])], floor(pars[5]), ceiling(pars[5]))
    
    
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3],pars[4],pars[5])
    
    # mixture weights for floor(K),ceiling(K) from real part of non-integer K
    realK <- parscont[3] - floor(parscont[3])
    fmW <- c(1-realK,realK)
    
  } else if (model == c("MK_P_RNplus")){

    
    K_range <- c(0:max(set_sizes))
    poissW <- dpois(c(0:(max(K_range)-1)), pars[5], log = FALSE)
    poissW <- c(poissW,1-sum(poissW))
    
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3],pars[4])
    
  } else if (model == "MK_U_RNplus") {
    
    K_range <- c(0:max(set_sizes))
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3],pars[4],pars[5])
    
  } else {
    
    K_range <- set_sizes
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3])
    
    if (grepl("RNplus",model)){
      parscont <- c(parscont,pars[4])
    }
    
    
  }
  
  
  out <- vector("numeric", length(set_sizes))
  
  for (i in seq_along(error_list)) {
    
    if (model %in% c("MK_F_RNplus")){
      
      out[[i]] <- numintroutineF(precision = precision, parscont = parscont, K_range = K_range, errors = error_list[[i]],
                                 set_sizes = set_sizes, sz = i)
      
    } else if (model %in% c("MK_FM_RNplus")){
      
      out[[i]] <- numintroutineFM(precision = precision, parscont = parscont, weights = fmW, K_range = K_range, errors = error_list[[i]],
                                 set_sizes = set_sizes, sz = i)
      
    } 
    
    else if (model %in% c("MK_P_RNplus")) {
      
      out[[i]] <- numintroutineP(precision = precision, parscont = parscont, weights = poissW, K_range = K_range, errors = error_list[[i]],
                                 set_sizes = set_sizes, sz = i)
      
    } else if (model %in% c("MK_U_RNplus")) {
      
      out[[i]] <- numintroutineU(precision = precision, parscont = parscont, K_range = K_range, errors = error_list[[i]],
                                 set_sizes = set_sizes, sz = i)
      
    } else {
      
      out[[i]] <- numintroutine(precision = precision[i], parscont = parscont, errors = error_list[[i]], model = model)
      
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

numintroutine <- function(precision, parscont, errors, model) {
  
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

numintroutineF <- function(precision, parscont, K_range, errors, set_sizes, sz) {
  
  coreFunction <- "cint_fun_MK_RNplus"
  
  pEncode <-  min(parscont[3],set_sizes[[sz]]) / set_sizes[[sz]]
  pars <- c(precision[[pmin(sz, length(precision))]],parscont[c(1,2)])
  
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
numintroutineFM <- function(precision, parscont, weights, K_range, errors, set_sizes, sz) {
  
  coreFunction <- "cint_fun_MK_RNplus"
  
  out <- vector("numeric", length(errors))
  pEncode <-  min(parscont[3],set_sizes[[sz]]) / set_sizes[[sz]]
  
  if (pEncode < 1){
    
    for (FMind in seq_along(weights)){
      
      err<- vector("numeric",length(errors))
      pars <- c(tail(precision,2)[[FMind]],parscont[c(1,2)])
      
      for (i in seq_along(err)) {
        err[i] <- vp_integration(error = errors[i], pars = pars, 
                                 coreFunction = coreFunction)
        
      }
      out <- out + weights[[FMind]] * (pEncode*err + (1-pEncode)*one_two_pi)
    }
    
    
  } else {
    
    
    err<- vector("numeric",length(errors))
    pars <- c(precision[[sz]],parscont[c(1,2)])
    
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

numintroutineU <- function(precision, parscont, K_range, errors, set_sizes, sz) {
 
  ## calculate proportion of real part of K from total K
  unifw_real <- (parscont[3] - floor(parscont[3])) / parscont[3] 
  ## calculate weight for integer K values from 0 to floor(K)
  unifw_int <- rep((1 - unifw_real)/ceiling(parscont[3]), ceiling(parscont[3]))
  ## weights for 0:ceiling(K)
  unifw_all <- c(unifw_int, unifw_real)
  ## restrict weights to current set_size and give remainder to last set_size
  unifw <- unifw_all[1:min(ceiling(parscont[3])+1,set_sizes[[sz]]+1)]
  unifw[length(unifw)] <- 1-sum(unifw[1:(length(unifw)-1)])

  coreFunction <- "cint_fun_MK_RNplus"
  
  out <- vector("numeric", length(errors))
  for (j in c(1:length(unifw))){ # K = 0 to K = min(floor(Kpar),setsizes[[sz]]) integer items
    
    K <- j - 1
    
    
    if (K < set_sizes[[sz]]){
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[j]],parscont[c(1,2)])
      
      for (i in seq_along(err)) {
        
        err[i] <- vp_integration(error = errors[i], pars = pars, 
                             coreFunction = coreFunction)
      }
      
      pEncode <- K/set_sizes[[sz]]
      out <- out +  unifw[[j]]*(pEncode*err + (1-pEncode)*one_two_pi)
      
    } else {
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[(set_sizes[[sz]] +1)]],parscont[c(1,2)])
      
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

numintroutineP <- function(precision, parscont, weights, K_range, errors, set_sizes, sz) {
  
  coreFunction <- "cint_fun_MK_RNplus"
  
  out <- vector("numeric", length(errors))
  for (j in c(1:length(weights))){ # K = 0 to K = max(setsize) integer items
    
    K <- j - 1 #to start from Kitems = 0
    
    
    if (K < set_sizes[[sz]]){
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[j]],parscont[c(1,2)])
      
      for (i in seq_along(err)) {
        
        err[i] <- vp_integration(error = errors[i], pars = pars, 
                                 coreFunction = coreFunction)
      }
      
      pEncode <- K/set_sizes[[sz]]
      out <- out +  weights[[j]]*(pEncode*err + (1-pEncode)*one_two_pi)
      
    } else {
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[(set_sizes[[sz]] + 1)]],parscont[c(1,2)])
      
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
