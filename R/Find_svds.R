.matprod <- function(v,args) args[[1]] %*% v
.tmatprod <- function(v,args) t(t(v) %*% args[[1]])
Find_svds <- function(A){
  RSpectra::svd(A=.matprod,k=ncol(A),nu=ncol(A),nv=0,
                 dim=dim(A),args=list(A),Atrans=.tmatprod)$u
}