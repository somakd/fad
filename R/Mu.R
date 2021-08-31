Mu <- function(x,Mubar,L,D0){
# Purpose: Compute unconstrained mu_mle as a function of x
  q = ncol(L)
  D = as.numeric(D0)
  mS = 1./(1+x*D)*Mubar - x*1/(1+x*D)*L%*%
    solve(diag(q) + x*crossprod(L, (1/(1+x*D)*L)))%*%
    crossprod(L, (1/(1+x*D)*Mubar))
  mu = sum(mS^2) - 1
  return(mu)
}
