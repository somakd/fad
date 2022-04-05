
AECM_cycle <- function(X, M, L, D0) {
# Purpose: One ECM cycle for ML method
# for N by P data matrix X
  
  N = nrow(X);                # number of observations
  P = nrow(L);                # number of coordinates
  q = ncol(L);                # number of factors
  
  if(is.null(q)) {
    q = 1
    L = as.matrix(L)
  }
  D = as.numeric(D0)
  
  # Expectation Step
  cmatdg = cmdg(L,D)

  val =  taom(X,L,D,M,cmatdg)
  Tao = val$tao
  Tao[Tao<1e-5] = 1e-5
  iSxM = val$iSxM
  m = iSxM/Tao
  sd = sqrt(1/Tao)
  ER = ERF(P,m,sd)
  ERR = P*sd^2 + m*ER
  
  # Conditional Maximization Step 
  
  # Mean -- Mu  
  
  Mubar = RXM(X,ER)
  if (sum(Mubar^2) >= 1){
      upper = sum(Mubar*1/D*Mubar - Mubar*1/D*L%*%(cmatdg*crossprod(L, 1/D*Mubar)))
      int0 = c(0,upper)
  }else{
    e = Find_Leig(L,D)
    b = abs(sum(e$vectors*Mubar))
    int0 = c((0.5*b-1)/e$values,0)
  }
  # x = fzero(Mu,Mubar=Mubar,L=L,D=D,x=int0)$x
  x = uniroot(f = Mu,interval = int0, Mubar=Mubar, L=L,D=D)
  
  mmu = 1./(1+x*D)*Mubar - x*1/(1+x*D)*L%*%
    solve(diag(q) + x*crossprod(L, (1/(1+x*D)*L)))%*%
    crossprod(L, (1/(1+x*D)*Mubar))
  mmu = c(mmu)
  #print(sum(mmu^2))
  
  # Update Expectations
  
  iSxM = ism(X,L,D,mmu,cmatdg)
  m = iSxM/Tao
  ER = ERF(P,m,sd)
  ERR = P*sd^2 + m*ER
  VR = ERR - ER^2
  VR[VR<0] = 0
  
  # Correlation parameters
  ybar = RXM(X,ER)
  csd = mSD(X,mmu,D,ER,VR,ybar)
  start = (1/csd^2)*D
  start[start >= 1] = 1
  zz = fads.fit.rxmd(X,mmu,csd,ER,VR,ybar,q, start = start, lower = 1e-5)
  
  return(list(converge = zz$converged, gerr = zz$gerr, loglik = zz$loglik,
              mmu = mmu, mvar = csd^2, mD = zz$uniquenesses, mL = zz$loadings))
}
