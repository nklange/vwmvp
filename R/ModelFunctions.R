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