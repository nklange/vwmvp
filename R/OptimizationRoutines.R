
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
  
  if (model == c("MK_F_RNplus")){
    
    K_range <- set_sizes[set_sizes <= pmin(max(set_sizes), pars[5])]
    # K_range <- K_range[K_range <= max(set_sizes)]
    # K_range <- K_range[K_range <= pars[5]]
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3],pars[4],pars[5])
    
  } else if (model == c("MK_P_RNplus")){
    
    K_range <- c(0:max(set_sizes))
    poissW <- dpois(c(0:(max(K_range)-1)), pars[5], log = FALSE)
    poissW <- c(poissW,1-sum(poissW))
    
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3],pars[4])
    
  }  else if (model == "MK_U_RNplus") {
    
    K_range <- c(0:max(set_sizes))
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3],pars[4],pars[5])
  }
  else {
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
      
    } else if (model %in% c("MK_P_RNplus")) {
      
      out[[i]] <- numintroutineP(precision = precision, parscont = parscont, poissW = poissW, K_range = K_range, errors = error_list[[i]],
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
  maxEnc <- min(parscont[3],set_sizes[[sz]])
  
  out <- vector("numeric", length(errors))
  
  err<- vector("numeric",length(errors))
  #pars <- c(precision[[match(maxEnc,K_range)]],parscont[c(1,2)])
  pars <- c(precision[[pmin(sz, length(precision))]],parscont[c(1,2)])
  
  for (i in seq_along(err)) {
    err[i] <- vp_integration(error = errors[i], pars = pars, 
                             coreFunction = coreFunction)
  }
  
  pEncode <- maxEnc/set_sizes[[sz]]
  out <- pEncode*err + (1-pEncode)*one_two_pi
  
  # out <- out/length(K_range)
  
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
  for (j in c(1:length(unifw))){ #runs from 0
    
    K <- j - 1 #to start from Kitems = 0
    
    
    if (K < set_sizes[[sz]]){
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[j]],parscont[c(1,2)])
      
      for (i in seq_along(err)) {
        
        err[i] <- vp_integration(error = errors[i], pars = pars, 
                             coreFunction = coreFunction)
      }
      
      pEncode <- K/set_sizes[[sz]]
      out <- out + unifw[[j]]*(pEncode*err + (1-pEncode)*one_two_pi)
      
    } else {
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[length(precision)],parscont[c(1,2)])
      
      for (i in seq_along(err)) {
        
        err[i] <- vp_integration(error = errors[i], pars = pars, 
                             coreFunction = coreFunction)
      }
      
      out <- out + unifw[[j]]*err
    }
    
    
  }
  # out <- out/length(K_range)
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
  
}

numintroutineP <- function(precision, parscont, poissW, errors, K_range, set_sizes, sz) {
  
  coreFunction <- "cint_fun_MK_RNplus"
  
  out <- vector("numeric", length(errors))
  for (j in c(1:length(K_range))){
    
    
    
    
    if (K_range[[j]] < set_sizes[[sz]]){
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[[j]],parscont)
      
      for (i in seq_along(err)) {
        
        err[i] <- tryCatch(
          integrate(match.fun(coreFunction), 
                    0, Inf, 
                    pars,
                    radian = errors[i], stop.on.error = FALSE)$value, 
          error = function(e) tryCatch(
            integrate(match.fun(coreFunction), 
                      0, Inf, 
                      pars,
                      radian = if (errors[i] == 0) {
                        circular(.Machine$double.xmin)  
                      } else errors[i], stop.on.error = FALSE)$value), 
          error = function(e) tryCatch(
            integrate(match.fun(coreFunction), 
                      0, Inf,   
                      pars,
                      radian = if (errors[i] == 0) {
                        circular(.Machine$double.eps^2) 
                      } else errors[i], stop.on.error = FALSE)$value), error = function(e) NA)
      }
      
      pEncode <- min(K_range[[j]]/set_sizes[[sz]],1)
      out <- out + poissW[[j]]*(pEncode*err + (1-pEncode)*(1/2/pi))
    } else {
      
      err<- vector("numeric",length(errors))
      pars <- c(precision[match(set_sizes[[sz]],K_range)],parscont)
      
      for (i in seq_along(err)) {
        
        err[i] <- tryCatch(
          integrate(match.fun(coreFunction), 
                    0, Inf, 
                    pars,
                    radian = errors[i], stop.on.error = FALSE)$value, 
          error = function(e) tryCatch(
            integrate(match.fun(coreFunction), 
                      0, Inf, 
                      pars,
                      radian = if (errors[i] == 0) {
                        circular(.Machine$double.xmin)  
                      } else errors[i], stop.on.error = FALSE)$value), 
          error = function(e) tryCatch(
            integrate(match.fun(coreFunction), 
                      0, Inf,   
                      pars,
                      radian = if (errors[i] == 0) {
                        circular(.Machine$double.eps^2) 
                      } else errors[i], stop.on.error = FALSE)$value), error = function(e) NA)
      }
      
      out <- out + poissW[[j]]*err
    }
    
    
  }
  # out <- out/length(K_range)
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
  
}
