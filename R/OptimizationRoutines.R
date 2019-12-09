
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
  
  precision <- pars[1]/(set_sizes^pars[2])
  parscont <- c(pars[3])
  
  if (grepl("RNplus",model)){
    parscont <- c(parscont,pars[4])
  }
  
  out <- vector("numeric", length(set_sizes))
  
  for (i in seq_along(error_list)) {
    out[[i]] <- numintroutine(pars = c(precision[i],parscont),errors=error_list[[i]],model=model)
  }
  
  return(sum(out))
}

