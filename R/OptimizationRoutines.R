
ll_vp_full <- function(pars, model, error_list, set_sizes, type ,...) {
  
  if (type == "sim") {
    ll_fun <- match.fun(paste0("csim_fun_",model))
  } else if (type == "numint") {
    ll_fun <- ll_vp_numint
  }
  
  precision <- pars[1]/(set_sizes^pars[2])
  parscont <- c(pars[3])
  
  if (grepl("RNplus",model)){
    parscont <- c(parscont,pars[4])
  }
  
  out <- vector("numeric", length(set_sizes))
  for (i in seq_along(error_list)) {
   
    out[[i]] <- ll_fun(pars = c(precision[i],parscont), 
                       errors = error_list[[i]],
                       model = model)
  }
  
  return(sum(out))
  
}

ll_vp_numint <- function(pars, errors, model) {


  out <- vector("numeric", length(errors))
  for (i in seq_along(out)) {
    out[i] <- tryCatch(
      integrate(match.fun(paste0("cint_fun_",model)), 
                0, Inf, 
                pars, 
                radian = errors[i], stop.on.error = FALSE)$value, 
      error = function(e) tryCatch(
        integrate(match.fun(paste0("cint_fun_",model)), 
                  0, Inf, 
                  pars,
                  radian = if (errors[i] == 0) {
                    circular(.Machine$double.xmin)  
                  } else errors[i], stop.on.error = FALSE)$value), 
      error = function(e) tryCatch(
        integrate(match.fun(paste0("cint_fun_",model)), 
                  0, Inf,   
                  pars,
                  radian = if (errors[i] == 0) {
                    circular(.Machine$double.eps^2) 
                  } else errors[i], stop.on.error = FALSE)$value), error = function(e) NA)
  }
  
  
  if (any(out == 0) | any(!is.finite(out))){
    
    return(1e6)
    
  } else {
    return(-sum(log(out)))
  }
  
}
