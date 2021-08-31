random.init <- function(inputs,q,init){
  # for random initialization of M,L,D, controlled by set.seed()

 # old.seed = .Random.seed
  set.seed(init)
  N = nrow(inputs)
  P = ncol(inputs)
  M = rnorm(P)
  M = 1/sqrt(sum(M^2))*M
  L = matrix(rnorm(P*q),P,q)
  D = runif(P,0.2,0.8)
  SD = sqrt(rowSums(L^2)+D)
  G = crossprod(L,1/D*L)
  Q = eigen(G)$vectors
  L = L%*%Q

 # set.seed(old.seed)
  return(list(M=M,L=L,D=D,SD=SD))
}
