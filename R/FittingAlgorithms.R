#' @import tidyverse
#' @import circular
#' @import GA


FitVP <- function(data, model, method, rep = 10, seqrun = 5, nsim = 1500, startpar = NULL) {
  
  if (method == "sim"){
    if (grepl("RNplus",model)){
      lower <- rep(.Machine$double.eps,4)
      upper <- c(400,5,150,200)
    } else if (grepl("RNminus",model)){
      lower <- rep(.Machine$double.eps,3)
      upper <- c(200,5,150)
    }
    
    res <- fit_one_vp_ga(data, rep, model, lower, upper, seqrun, nsim, parallel = FALSE)
    
  } else if (method == "numint"){
    
    res <- fit_one_vp_nlminb(data, rep, model, startpar = NULL)
  }
  return(res)
}



# Genetic Algorithm ----------------------------------------------------------------------------------
fit_one_vp_ga <- function(data, rep, model, lower, upper, parallel = parallel, nsim = nsim, seqrun = seqrun,...) {
  
  dp <- prep_data_index(data)
  
  out_list <- vector("list", rep)
  
  for (i in seq_len(rep)) {
    tic <- Sys.time()
    tmp <- tryCatch(ga(type = "real-valued", fitness = ll_vp_sim,
                       error_list = dp$datalist, set_sizes = dp$set_sizes,
                       model = model, nsim = nsim,
                       ...,
                       #monitor = FALSE,
                       lower = lower, upper = upper, 
                       run = seqrun, parallel = parallel), 
                    error = function(e) NA)
    if (inherits(tmp, "ga")) {
      pars <- as_tibble(tmp@solution)
      colnames(pars) <- names(get_start_vp(model))
      pars$objective <- tmp@fitnessValue
      pars$iter <- tmp@iter
      pars$model <- model
      pars$time <- Sys.time() - tic
      pars$id <- dp$id
      pars$exp <- dp$exp
      pars$leftout <- dp$leftout
      pars$rep <- i
      pars$cvid <- dp$cvid
      out_list[[i]] <- pars
    }
  }
  bind_rows(out_list)
}


ll_vp_sim <- function(pars, model, error_list, set_sizes, nsim...){
  
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


# Numerical Integration NLMINB ----------------------------------------------------------

fit_one_vp_nlminb <- function(model, data, startpar, rep, ...) {
  
  dp <- prep_data(data)
  out_list <- vector("list", rep)
  
  if (is.null(startpar)){
    
  } else {
  startpar <- data.frame(do.call(rbind, startpar))
  startpar <- startpar %>% 
    filter(cvid == dp$cvid) %>% 
    select(names(get_start_vp(model))) 
  
  }
  
  for (i in seq_len(rep)) {
    
    start <- startpar %>% 
      slice(i) %>% 
      unlist()

    tic <- Sys.time()
    
    tmp <- tryCatch(nlminb(start, objective = ll_vp_numint, 
                           error_list = dp$datalist, set_sizes = dp$set_sizes, 
                           model = model,..., 
                           lower = .Machine$double.eps, upper = Inf,
                           control = list(
                             eval.max = 300, 
                             iter.max = 300, 
                             trace = 1,
                             rel.tol = 1e-5, ## default 1e-10, at 1e-3 difference between runs increases too much
                             x.tol = 1.5e-5 ## default 1.5e-8
                           )), error = function(e) NA)
    if (is.list(tmp)) {
      
      out_list[[i]] <- bind_cols(
        spread(enframe(tmp$par), name, value) %>%
          setNames(sort(names(get_start_vp(model))))
        , as_tibble(tmp[c(2, 3, 4, 6)])
        , tibble(
          time = Sys.time() - tic
        )
        , tibble (
          model = model
        )
        , tibble (
          id = dp$id
        )
        ,tibble (
          leftout = dp$leftout
        )
        , tibble (
          rep = i
        )
        ,tibble (
          exp = dp$exp
        )
        , tibble (
          cvid = dp$cvid
        )
        , spread(enframe(start, "start"), start, value) %>%
          setNames(paste0("s_",sort(names(get_start_vp(model)))))
      )
    } 
  }
  bind_rows(out_list)
}
