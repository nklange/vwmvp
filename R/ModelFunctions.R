sim_fun_MK_RNminus<- function(x,pars,base_radians){
  
  kappas <- rgamma(x, shape = pars[1] / pars[2], scale = pars[2])
  
  pchoose <- vector("numeric", 181)
  
  for (i in seq_len(x)) {
    
    pchoose <- pchoose +
      1 / (2 * pi * besselI(
        x = kappas[i],
        nu = 0,
        expon.scaled = TRUE
      )) *
      (exp(cos(base_radians) - 1)) ^ kappas[i]
    
  }
  
  return((pchoose / sum(pchoose)) * 0.5)
}


sim_fun_MK_RNplus<- function(x, pars, base_radians){
  
  kappa_r <- pars[3]
  
  
  kappas <- rgamma(x, shape = pars[1] / pars[2], scale = pars[2])
  
  pchoose <- vector("numeric", 181)
  kc <- vector("numeric",181)
  
  for (i in seq_len(x)) {
    
    kc <- sqrt(kappa_r^2 + kappas[i]^2 + 2 * kappa_r*kappas[i]*cos(base_radians))
    
    pchoose <- pchoose + 
      ((besselI(kc,0,expon.scaled = TRUE) / 
          (2*pi*besselI(kappas[i],0,expon.scaled = TRUE) * 
             besselI(kappa_r,0,expon.scaled = TRUE))) *
         exp(kc - (kappas[i] + kappa_r)))
    
  }
  
 return((pchoose / sum(pchoose)) * 0.5)
}


sim_fun_J_RNminus <- function(x, pars,base_radians){
  
  
  J <- rgamma(x, shape = pars[1]/pars[2], scale = pars[2])
  
  kappas <- KappafromJ(J)
  
  pchoose <- vector("numeric", 181)
  
  for (i in seq_len(x)) {
    
    pchoose <- pchoose +
      1/(2 * pi * besselI(x = kappas[i], nu = 0, expon.scaled = TRUE)) * 
      (exp(cos(base_radians) -1))^kappas[i]
  }
  
  return((pchoose/sum(pchoose) ) * 0.5)
}

sim_fun_J_RNplus <- function(x, pars, base_radians){
  
  kappa_r <- pars[3]
  J <- rgamma(x, shape = pars[1]/pars[2], scale = pars[2])
  kappas <- KappafromJ(J)
  
  pchoose <- vector("numeric", 181)
  kc <- vector("numeric",181)
  
  for (i in seq_len(x)) {
    
    kc <- sqrt(kappa_r^2 + kappas[i]^2 + 2 * kappa_r*kappas[i]*cos(base_radians))
    
    pchoose <- pchoose + 
      ((besselI(kc,0,expon.scaled = TRUE) / 
          (2*pi*besselI(kappas[i],0,expon.scaled = TRUE) * 
             besselI(kappa_r,0,expon.scaled = TRUE))) *
         exp(kc - (kappas[i] + kappa_r)))
    
  }
  
  return((pchoose/sum(pchoose)) * 0.5)
}
