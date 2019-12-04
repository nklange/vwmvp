#' @import tidyverse
#' @import circular
#' @import GA


FitVP <- function(data, model, type, rep = 10, startpar = NULL) {
  if (type == "sim"){
    if (grepl("RNplus",model)){
      lower <- lowerRNplus
      upper <- upperRNplus
    } else if (grepl("RNminus",model)){
      lower <- lowerRNminus
      upper <- upperRNminus
    }
    res <- fit_one_vp_ga(data, rep, model, type, lower, upper, parallel = FALSE)
  } else if (type == "numint"){
    res <- fit_one_vp_nlminb(data, rep, model, type, startpar = NULL)
  }
  return(res)
}



# Genetic Algorithm ----------------------------------------------------------------------------------


fit_one_vp_ga <- function(data, rep, model, type, lower, upper, ..., parallel = FALSE) {
  dp <- prep_data_index(data)
  
  out_list <- vector("list", rep)
  
  for (i in seq_len(rep)) {
    tic <- Sys.time()
    tmp <- tryCatch(ga(type = "real-valued", fitness = ll_vp_full,
                       error_list = dp$datalist, set_sizes = dp$set_sizes,
                       model = model,
                       ...,
                       #monitor = FALSE,
                       lower = lower, upper = upper, 
                       run = 5, parallel = parallel), 
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
    
    tmp <- tryCatch(nlminb(start, objective = ll_vp_full, 
                           type = type,
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
