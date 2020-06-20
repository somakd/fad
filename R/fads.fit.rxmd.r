
fads.fit.rxmd <-
  function(X,mmu,csd,ER,VR,ybar,q,
           start=NULL, maxit = 500, iSD=NULL, mu = NULL,lower = 0.005, control = NULL, ...)  {
    p <- length(mmu)
    if(is.null(start))
    {
      eigres <- eigs_sym_RXmD(X,mmu,1/csd,ER,VR,ybar,q,1)
      start <- 1 - rowSums( .postmdiag(eigres$vectors,sqrt(eigres$values))^2);
      start[start <= 0] = lower;
      start = cbind(start)
    }
    # print(start)
    fngr.fads <- function(Psi) FAfngrC.fads(Psi,X,mmu,csd,ER,VR,ybar,q)

    res <- optim_fad(par = start,fngr = fngr.fads,lower=lower,upper=1,method="L-BFGS-B",
                       control = c(list(fnscale=1,parscale = rep(0.01, length(start)),maxit = maxit,
                                        factr=1e2)))

    Lambda <- FAoutC.fads(res$par, X,mmu,csd,ER,VR,ybar, q,iSD,mu)
    dof <- 0.5 * ((p - q)^2 - p - q)
    #un <- setNames(res$par, colnames(R))
    class(Lambda) <- "loadings"
    gerr = sqrt(sum({rowSums(Lambda^2)+res$par - 1}^2))
    ans <- list(converged = res$convergence == 0,
                loadings = Lambda, uniquenesses = res$par,
                gerr = gerr,
                criteria = c(objective = res$value, counts = res$counts),
                factors = q, dof = dof, method = "mle")
    class(ans) <- "fad"
    return(ans);
  }


FAfngrC.fads <- function(Psi,X,mmu,csd,ER,VR,ybar, q)
{
  sc <- 1/sqrt(Psi)
  eigsres <- eigs_sym_RXmD(X,mmu,sc/csd,ER,VR,ybar,q,1)
  L <- eigsres$vectors

  e <- eigsres$values
  e <-  -sum(log(Psi)) - sum(1/Psi) -sum(log(e) - e) - q

  load <- .postmdiag(L , sqrt(pmax(eigsres$values-1,0)))
  load <- {1/sc}*load

  g <- rowSums(load^2) + Psi - 1;
  g <- g/Psi^2;
  return(list(-e,g));
}


FAoutC.fads <- function(Psi, X,mmu,csd,ER,VR,ybar, q,iSD,mu)
{
  sc <- 1/sqrt(Psi)
  eigsres <- eigs_sym_RXmD(X,mmu,sc/csd,ER,VR,ybar,q,1)
  load <- .postmdiag(eigsres$vectors , sqrt(pmax(eigsres$values-1,0)))
  load <- {1/sc}*load

  return(load)
}


