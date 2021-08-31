
cppFunction('NumericMatrix postmdiag(NumericMatrix X,NumericVector d) {
  int c = X.ncol();
  int r = X.nrow();
  double dj;
  if(d.size() != c)
  {
    Rcpp::stop("Length of d must be same as number of columns of X"); 
  }
  for(int j=0;j<c;++j)
  {
    dj = d[j];
    for(int i=0;i<r;++i)
    X[j*r+i] *= dj;
  }
  return(X);
}')