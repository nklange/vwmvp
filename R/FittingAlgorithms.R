

FitVP <- function(data, model, method, rep = rep, seqrun = 5, nsim = 1500, startpar = NULL) {
  
  
  res <- NULL
  
  if (method == "sim"){
    
    lower <- NULL
    upper <- NULL
    if (grepl("RNplus",model)){
      lower <- rep(.Machine$double.eps,4)
      upper <- c(400,5,150,200)
    } else if (grepl("RNminus",model)){
      lower <- rep(.Machine$double.eps,3)
      upper <- c(200,5,150)
    }
    
    res <- fit_one_vp_ga(data = data, rep = rep, model = model, lower = lower, upper = upper, 
                         nsim = nsim, seqrun = seqrun)
    
    
    
  } else if (method == "numint"){
    
    res <- fit_one_vp_nlminb(data = data, rep = rep, model = model, startpar = startpar)
    
  }
  
  return(res)
}

# Genetic Algorithm ----------------------------------------------------------------------------------

# Function

fit_one_vp_ga <- function(data = data, rep = rep, model = model, lower = lower, upper = upper, 
                          nsim = nsim, seqrun = seqrun, parallel = FALSE) {
  dp <- prep_data_index(data)
  
  out_list <- vector("list", rep)
  
  for (i in seq_len(rep)) {
    tic <- Sys.time()
    tmp <- tryCatch(GA::ga(type = "real-valued", fitness = ll_vp_sim,
                       error_list = dp$datalist, set_sizes = dp$set_sizes,
                       model = model,
                       ...,
                       #monitor = FALSE,
                       lower = lower, upper = upper, nsim = nsim,
                       run = seqrun, parallel = parallel), 
                    error = function(e) NA)
    if (inherits(tmp, "ga")) {
      pars <- tibble::as_tibble(tmp@solution)
      colnames(pars) <- names(get_start_vp(model))
      pars$objective <- tmp@fitnessValue
      pars$iter <- tmp@iter
      pars$type <- tmp@type
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
  dplyr::bind_rows(out_list)
}

# Numerical Integration NLMINB ----------------------------------------------------------
fit_one_vp_nlminb <- function(data, rep, model, startpar, ...) {
  
  dp <- prep_data(data)
  out_list <- vector("list", rep)
  
  if (is.null(startpar)){
    startpar <- NULL
    for (reps in c(1:rep)){
      startpar[[reps]] <- get_start_vp(model)
    }
    startpar <- data.frame(matrix(unlist(startpar), nrow=length(startpar), byrow=T))
    colnames(startpar) <- names(get_start_vp(model))
    
  } else {
    startpar <- data.frame(do.call(rbind, startpar))
    startpar <- startpar %>% 
      filter(cvid == dp$cvid) %>%
      select(names(get_start_vp(model))) 
    
  }
  
  for (i in seq_len(rep)) {
    #print(startpar)
    
    start <- startpar %>% 
      slice(i) %>% 
      unlist()
    
    
    tic <- Sys.time()
    
    tmp <- tryCatch(nlminb(start, objective = ll_vp_numint, 
                           error_list = dp$datalist, set_sizes = dp$set_sizes, 
                           model = model, ..., 
                           lower = .Machine$double.eps, upper = Inf,
                           control = list(
                             eval.max = 300, 
                             iter.max = 300, 
                             trace = 1,
                             rel.tol = 1e-5, ## default 1e-10, at 1e-3 difference between runs increases too much
                             x.tol = 1.5e-8 ## default 1.5e-8
                           )), error = function(e) NA)
    if (is.list(tmp)) {
      
      out_list[[i]] <- bind_cols(
        tidyr::spread(tibble::enframe(tmp$par), name, value) %>%
          setNames(sort(names(get_start_vp(model))))
        , tibble::as_tibble(tmp[c(2, 3, 4, 6)])
        , tibble::tibble(
          time = Sys.time() - tic
        )
        , tibble::tibble (
          model = model
        )
        , tibble::tibble (
          id = dp$id
        )
        , tibble::tibble (
          leftout = dp$leftout
        )
        , tibble::tibble (
          rep = i
        )
        , tibble::tibble (
          exp = dp$exp
        )
        , tibble::tibble (
          cvid = dp$cvid
        )
        , tidyr::spread(tibble::enframe(start, "start"), start, value) %>%
          setNames(paste0("s_",sort(names(get_start_vp(model)))))
      )
    } 
  }
  dplyr::bind_rows(out_list)
}