#' Factor Analysis for data on a sphere (high or low dimensional).
#' @rdname fads
#' @description Perform fast matrix-free maximum-likelihood factor analysis on
#'   data on sphere, works if number of variables is more than number of 
#'   observations.
#' @param inputs A numeric matrix or an object that can be coerced to a
#'   numeric matrix.
#' @param q The number of factors to be fitted.
#' @param ii The random seeds for initialization. Default 123 if no initial 
#'   values of parameters are imported.
#' @param M The initial values of mean.
#' @param L The initial values of loading matrix.
#' @param D The initial values of uniquenesses.
#' @param gamma The common constant used in the eBIC formula. Default 'NA'.
#' @param maxiter The maximum iterations. Default 10,000
#' @param epsi The absolute difference between final data log-likelihood values
#'   on consecutive step. Default 0.0001.
#'
#' @return An object of class \code{"fads"} with components
#' \item{mu}{The estimate mean.}
#'\item{loadings}{A matrix of loadings on the correlation scale, one column for each factor.  The
#' factors are ordered in decreasing order of sums of squares of
#' loadings, and given the sign that will make the sum of the loadings
#' positive.This is of class \code{"loadings"}}
#' \item{uniquenesses}{The uniquenesses computed on the correlation scale.}
#' \item{sd}{The estimated standard deviations.}
#' \item{iter}{The number of iterations}
#' \item{gerr}{the difference between the gradients on consecutive step.}
#' \item{loglik,eBIC}{The maximum log-likelihood the extended Bayesian 
#' Information Criteria (Chen and Chen,2008).}
#'
#' @export
fads = function(inputs,q,ii = 123,M=NULL,L=NULL,D=NULL,gamma=NA,
                                 maxiter = 10000,epsi = 1e-4){
  # Purpose: AECM FOR FACTOR ANALYSIS OF DIRECTIONAL DATA
  # output: mL,mD: correlation parameters
  N = nrow(inputs)
  P = ncol(inputs)            
  iter = 0                    # loop step

  # Initialize Mean & Covariance Components
  if(is.null(M)){
    init = random.init(inputs,q,ii)
    M = init$M
    L = init$L
    D = init$D
  }
  
  G = crossprod(L,1/D*L)
  Q = eigen(G)$vectors
  L = L%*%Q

  llk0 = LLk_Dir(inputs,M,L,D)
  
  print(c(iter,llk0))
  
  for (iter in 0:maxiter){
    par0 = c(CartToSph(M),L,log(D)) # \theta_0
    out = AECM_cycle(inputs,M,L,D)
    #print(out$converge)
    M = out$mmu
    L = sqrt(out$mvar)*out$mL
    D = out$mvar*out$mD
    iter = iter+1
    print(iter)
    par1 = c(CartToSph(M),L,log(D)) # \theta_1
    r = par1 - par0 # r 
    
    out = AECM_cycle(inputs,M,L,D)
    M = out$mmu
    L = sqrt(out$mvar)*out$mL
    D = out$mvar*out$mD
    iter = iter+1
    par2 = c(CartToSph(M),L,log(D)) # \theta_2
    v = par2 - par1 - r; # v
    
    alpha = -sqrt(sum(r^2))/sqrt(sum(v^2)) # steplength, Varadhan & Roland (2008), eq.9
    
    # Modify alpha, Varadhan & Roland (2008), Section 6
    if (alpha > -1) alpha = -1
    
    pnew = par0 - 2*alpha*r + alpha^2*v
    aM = SphToCart(pnew[1:(P-1)])
    aL = matrix(pnew[P:(P-1+P*q)],P,q)
    aD = exp(pnew[(P+P*q):length(pnew)])
    
    aG = crossprod(aL,1/aD*aL)
    aQ = eigen(aG)$vectors
    aL = aL%*%aQ
    
    llka = LLk_Dir(inputs,aM,aL,aD)
  
    while (llka < llk0){
      alpha = (alpha-1)/2
      if (alpha == -1){
        break
      }
      pnew = par0 - 2*alpha*r + alpha^2*v
      aM = SphToCart(pnew[1:(P-1)])
      aL = matrix(pnew[P:(P-1+P*q)],P,q)
      aD = exp(pnew[(P+P*q):length(pnew)])
      
      aG = crossprod(aL,1/aD*aL)
      aQ = eigen(aG)$vectors
      aL = aL%*%aQ
      
      llka = LLk_Dir(inputs,aM,aL,aD)
    }
   
    out = AECM_cycle(inputs,aM,aL,aD)
    #print(out$converge)
    M = out$mmu
    L = sqrt(out$mvar)*out$mL
    D = out$mvar*out$mD
    iter = iter+1
    llk = LLk_Dir(inputs,M,L,D)

    print(c(iter,llk0, llk, range(out$mvar),range(out$mD),range(out$mL)))

    #print(c(alpha,llk0,llk))
    if (abs(llk-llk0) < epsi){
      break
      }
    llk0 = llk

  }
  
  #GLS to estimate factors (scaled out by out$msd)
  #svd.H = svd(crossprod(out$mL,as.numeric(1/out$mD)*out$mL))
  #Z <- out$Y %*% (as.numeric(1/out$mD)*out$mL) %*% (svd.H$u %*% (1/svd.H$d * t(svd.H$v)))
  
  if(is.na(gamma)){
    gamma = max(1-1/(2*log(P,base = N)),0)
  }
  eBIC = -2*llk + P*q*(log(N) + 2*gamma*log(P))
  
  fit <- list(iter = iter, gerr = out$gerr, loglik = llk, eBIC = eBIC,
              mu = out$mmu, sd = sqrt(out$mvar), uniquenesses = out$mD, loadings = out$mL)
  class(fit$loadings) <- "loadings"
  class(fit) <- "fads"
  return(fit)
}
