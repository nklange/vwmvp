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


# preprocess data -----------------------------------------------------------

# make ppl get their data into right format?? ie both degrees and radians

# get parameter names/starting paras ----------------------------------------
# Get start parameters ---------------------------------------------

get_start_vp <- function(model) {
  
  
  
  #if (is.null(res_ga)) {
  if (model == "VP(K)A-") {
    start <- c(
      mkappa1 = runif(1, 30, 60),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 10, 40)
    )
  } else if (model == "VP(K)A+")  {
    start <- c(
      mkappa1 = runif(1, 100, 200),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 50, 130),
      kappa_r = runif(1, 10, 70)
    )
  } else if (model == "VP(J)A-") {
    start <- c(
      J1bar = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 10, 40)
    )
  } else if (model == "VP(J)A+") {
    start <- c(
      J1bar = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 10, 40),
      kappa_r = runif(1, 30, 60)
    )
    #  }
  } else if (model == "EP(K)A+") {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      kappa_r = runif(1, 30, 60)
    )
    #  }
  } else if (model %in% c("EP(K)F+","EP(K)P+","EP(K)U+","EP(K)FM+")) {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      kappa_r = runif(1, 30, 60),
      K = runif(1,0,10)
    )
  }  else if (model %in% c("EP(K)F-","EP(K)P-","EP(K)U-","EP(K)FM-")) {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      K = runif(1,0,10)
    )
  } else if (model %in% c("EP(K)A-")) {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2)
    )
  } else if (model %in% c("SA(K)F+","SA(K)P+","SA(K)U+")) {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      kappa_r = runif(1, 30, 60),
      K = runif(1,0,10)
    )
    #  }
  } else if (model %in% c("SA(K)F-","SA(K)P-","SA(K)U-")) {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      K = runif(1,0,10)
    )
    #  }
  } else if (model == "SA(K)A+") {
    start <- c(
      kappa_r = runif(1, 1, 30)
    )
    #  }
  } else if (model %in%  c("VP(K)P+","VP(K)U+",
                           "VP(K)F+","VP(K)FM+",
                           "VP(K)FM+","VP(K)P2+","VP(K)U2+"))  {
    start <- c(
      mkappa1 = runif(1, 100, 200),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 50, 130),
      kappa_r = runif(1, 10, 70),
      K = runif(1,0,10)
    )
  } else if (model %in%  c("VP(K)P-","VP(K)U-",
                           "VP(K)F-","VP(K)FM-",
                           "VP(K)FM-","VP(K)P2-","VP(K)U2-"))  {
    start <- c(
      mkappa1 = runif(1, 50, 100),
      alpha = runif(1, 0.5, 2),
      tau = runif(1, 50, 130),
      K = runif(1,0,10)
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