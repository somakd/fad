.LLDmatprod <- function(v,args) args[[1]] %*% crossprod(args[[1]],v) + args[[2]]*v
Find_Leig = function(L,D){
  eigs_sym(.LLDmatprod,1,n=nrow(L),args=list(L,D))
}