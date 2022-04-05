LLk_Dir <- function(X,mu0,L,D0){
# Purpose: Calculate Marginal Loglikehood Kernel for Directional Data
  mu = as.numeric(mu0)
  D = as.numeric(D0)
  N = nrow(X);            # number of obs.
  P = nrow(L);		  # dimensionality of the directional input
  q = ncol(L);            # number of factors
 
  cmatdg = cmdg(L,D)
  logdetS = sum(log(D)) + 2*sum(log(1/sqrt(cmatdg)))
  iSm = mu/D - {1/D}*L%*%(cmatdg*crossprod(L, mu/D))
  val =  taom(X,L,D,mu,cmatdg)
  Tao = val$tao
  Tao[Tao<1e-5] = 1e-5
  iSxM = val$iSxM
  m = iSxM/Tao
  sd = sqrt(1/Tao)
  
  logiv = logIv(P-1,m,sd);
  llk = - N*0.5*logdetS - N*0.5*sum(mu*iSm) + sum(0.5*Tao*m^2 + logiv);
  return(llk)
}
