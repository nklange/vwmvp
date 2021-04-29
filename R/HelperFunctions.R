# constants 

circ_null <- circular::circular(0)
base_radians <- circular::circular((0:180) * 2 * pi/360)

kappa_map <- c(seq(0, 10, length.out = 250),seq(10.001,1e4,length.out = 250))
J_map <- kappa_map*besselI(kappa_map,1,expon.scaled = TRUE)/besselI(kappa_map,0,expon.scaled = TRUE)
mapJkappa <- approxfun(J_map,kappa_map, yleft = 0)

KappafromJ <- function(J) {
  out <- mapJkappa(J)
  ifelse(is.na(out), J, out)
}

# Pewsey (2004)
circkurtosis <- function(errors){ 
  sum(cos(2*(errors - mean(errors))))/length(errors)
}

circkurtosis_std <-function(errors){
  rbar <-  sum(cos((errors - mean(errors))))/length(errors)
  alpha2 <- sum(cos(2*(errors - mean(errors))))/length(errors)
  (alpha2 - rbar^4) / ((1 - rbar)^2)
}

# preprocess data -----------------------------------------------------------

# make ppl get their data into right format?? ie both degrees and radians

# get parameter names/starting paras ----------------------------------------
# Get start parameters ---------------------------------------------

get_start_vp <- function(model) {
  
  
  
  #if (is.null(res_ga)) {
  if (model == "MK_RNminus") {
    start <- c(
      mkappa1 = runif(1, 30, 60),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 10, 40)
    )
  } else if (model == "MK_RNplus")  {
    start <- c(
      mkappa1 = runif(1, 100, 200),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 50, 130),
      kappa_r = runif(1, 10, 70)
    )
  } else if (model == "J_RNminus") {
    start <- c(
      J1bar = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 10, 40)
    )
  } else if (model == "J_RNplus") {
    start <- c(
      J1bar = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 10, 40),
      kappa_r = runif(1, 30, 60)
    )
    #  }
  } else if (model == "EP_RNplus") {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      kappa_r = runif(1, 30, 60)
    )
    #  }
  } else if (model %in% c("EP_F_RNplus","EP_P_RNplus","EP_U_RNplus","EP_FM_RNplus",
                          "EP_P2_RNplus","EP_U2_RNplus","EP_FM2_RNplus")) {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      kappa_r = runif(1, 30, 60),
      K = runif(1,0,10)
    )
  }  else if (model %in% c("EP_F_RNminus","EP_P_RNminus","EP_U_RNminus","EP_FM_RNminus",
                           "EP_P2_RNminus","EP_U2_RNminus","EP_FM2_RNminus")) {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      K = runif(1,0,10)
    )
  } else if (model %in% c("EP_RNminus")) {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2)
    )
  } else if (model %in% c("SA_F_RNplus","SA_P_RNplus","SA_U_RNplus")) {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      kappa_r = runif(1, 30, 60),
      K = runif(1,0,10)
    )
    #  }
  } else if (model %in% c("SA_F_RNminus","SA_P_RNminus","SA_U_RNminus")) {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      K = runif(1,0,10)
    )
    #  }
  } else if (model %in% c("SA_RNplus","VMnosetsize")) {
    start <- c(
      kappa_r = runif(1, 1, 30)
    )
    #  }
  } else if (model %in%  c("MK_P_RNplus","MK_U_RNplus","MK_F_RNplus","MK_FM_RNplus",
                           "MK_FM2_RNplus","MK_P2_RNplus","MK_U2_RNplus"))  {
    start <- c(
      mkappa1 = runif(1, 100, 200),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 50, 130),
      kappa_r = runif(1, 10, 70),
      K = runif(1,0,10)
    )
  } else if (model %in%  c("MK_P_RNminus","MK_U_RNminus","MK_F_RNminus","MK_FM_RNminus","MK_FM2_RNminus",
                           "MK_P2_RNminus","MK_U2_RNminus"))  {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 50, 130),
      K = runif(1,0,10)
    )
  } else if (model %in%  c("UVM"))  {
    start <- c(
      kuvm = runif(1, 0, 10),
      wuvm = runif(1,0, 1)
    )
  }  else if (model == "VPnosetsize") {
    start <- c(
      kappa = runif(1, 1, 50),
      tau = runif(1, 10, 40)
    )
  } else if (model == "VPplusnosetsize") {
    start <- c(
      kappa = runif(1, 100, 200),
      tau = runif(1, 10, 40),
      kappa_r = runif(1,10,70)
    )
  }
    
  return(start)
}

# nlminb: fit radians -------------------------------------------------------

prep_data <- function(data) {
  
  dl <- split(circular::circular(data$error_0), f = data$set_size)
  id <- unique(data$id)
  set_sizes <- sort(unique(data$set_size))
  
  if (is.null(unique(data$cvid))) {
    cvid <- id
  } else {
    cvid <- unique(data$cvid)
  }
  
  if (is.null(unique(data$exp))){
    exp <- NA
  } else {
    exp <- unique(data$exp)
  }
  
  if (is.null(unique(data$leftout))){
    leftout <- 0
  } else {
    leftout <- unique(data$leftout)
  }
  
  stopifnot(names(dl) == as.character(set_sizes))
  return(list(
    datalist = dl,
    set_sizes = set_sizes,
    id = id,
    cvid = cvid,
    exp = exp,
    leftout = leftout
  ))
}

# GA data: fitting degrees -------------------------------------------------
prep_data_index <- function(data) {
  
  dl <- split((abs(data$deg_error_0) + 1), f = data$set_size)
  id <- unique(data$id)
  set_sizes <- sort(unique(data$set_size))
  
  if (is.null(unique(data$cvid))) {
    cvid <- id
  } else {
    cvid <- unique(data$cvid)
  }
  
  if (is.null(unique(data$exp))){
    exp <- NA
  } else {
    exp <- unique(data$exp)
  }
  
  if (is.null(unique(data$leftout))){
    leftout <- 0
  } else {
    leftout <- unique(data$leftout)
  }
  stopifnot(names(dl) == as.character(set_sizes))
  return(list(
    datalist = dl,
    set_sizes = set_sizes,
    id = id,
    cvid = cvid,
    exp = exp,
    leftout = leftout
  ))
}

# Adjust precision and weight parameters -----------------------------------

prep_parameters <- function(pars, model, set_sizes){
  
  precision <- NULL
  parscont <- NULL
  weights <- NULL
  
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
    weights <- c(1-realK,realK)
    
  } else if (model %in% c("MK_P_RNplus","MK_P_RNminus","MK_P2_RNplus","MK_P2_RNminus")){
    
    
    K_range <- c(0:max(set_sizes))
    poissW <- dpois(c(0:(max(K_range)-1)), pars[length(pars)], log = FALSE)
    weights <- c(poissW,1-sum(poissW))
    
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3:length(pars)])
    
  } else if (model %in% c("MK_U_RNplus","MK_U_RNminus","MK_U2_RNplus","MK_U2_RNminus")) {
    
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
    
  } else if (model %in% c("EP_FM_RNplus","EP_FM_RNminus","EP_FM2_RNplus","EP_FM2_RNminus")){
    
    K_range <- c(set_sizes[set_sizes <= pmin(max(set_sizes), pars[length(pars)])], 
                 floor(pars[length(pars)]), ceiling(pars[length(pars)]))
    
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3:length(pars)])
    
    # mixture weights for floor(K),ceiling(K) from real part of non-integer K
    realK <- pars[length(pars)] - floor(pars[length(pars)])
    weights <- c(1-realK,realK)
    
    
  } else if (model %in% c("EP_P_RNplus","EP_P_RNminus","EP_P2_RNplus","EP_P2_RNminus")){
    
    
    K_range <- c(0:max(set_sizes))
    poissW <- dpois(c(0:(max(K_range)-1)), pars[length(pars)], log = FALSE)
    weights <- c(poissW,1-sum(poissW))
    
    precision <- pars[1]/(K_range^pars[2])
    
    if (grepl("RNplus",model)) {
      parscont <- c(pars[3:length(pars)])
    } else {
      parscont <- c(pars[length(pars)])
    }
    
    
  } else if (model %in% c("EP_U_RNplus","EP_U_RNminus","EP_U2_RNplus","EP_U2_RNminus")) {
    
    K_range <- c(0:max(set_sizes))
    precision <- pars[1]/(K_range^pars[2])
    parscont <- c(pars[3:length(pars)])
    
  } else if (model %in% c("SA_RNplus","VMnosetsize")){
    
    parscont <- pars[1]
    
  } else if (model %in% c("SA_F_RNplus","SA_U_RNplus","SA_F_RNminus","SA_U_RNminus","VPnosetsize")){
    
    precision <- pars[1]
    parscont <- c(pars[2:length(pars)])
    
  }  else if (model %in% c("SA_P_RNplus","SA_P_RNminus")){
    
    K_range <- c(0:max(set_sizes))
    poissW <- dpois(c(0:(max(K_range)-1)), pars[length(pars)], log = FALSE)
    weights <- c(poissW,1-sum(poissW))
    
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
    
    
  } else if (model == "UVM"){
    parscont <- pars
  }
  
  else if (model == "VPnosetsize"){
    precision <- pars[1]
    parscont <- pars[2]
    
  }
  else if (model == "VPplusnosetsize"){
    precision <- pars[1]
    parscont <- c(pars[2],pars[3])
    
  }
  return(list(precision,parscont,weights))
  
}