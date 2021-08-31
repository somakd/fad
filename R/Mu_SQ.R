Mu_SQ <- function(x,Mubar,D0){
  # Purpose: Compute unconstraited mu_mle as a function of x
  D = as.numeric(D0)
  
  mS = 1./(1+x*D)*Mubar 
  mu = sum(mS^2) - 1
  return(mu)
}
