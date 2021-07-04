
#ifndef _UTILS_CC_
#define _UTILS_CC_

#define STRICT_R_HEADERS
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export(.postmdiag)]]
NumericMatrix postmdiag(NumericMatrix X,NumericVector d)
{
  int c = X.ncol();
  int r = X.nrow();
  double dj;
  NumericMatrix Y(r,c);

  if(d.size() != c)
  {
    Rcpp::stop("Length of d must be same as number of columns of X");
  }
  for(int j=0;j<c;++j)
  {
    dj = d[j];
    NumericMatrix::iterator iX = X.begin()  + j*r;
    NumericMatrix::iterator iY = Y.begin()  + j*r;
    for(int i=0;i<r;++i)
      *(iY++) = *(iX++)*dj;
  }
  return(Y);
}


/*
 *The followuing function computes the 1/rowSum( (D*(X- 1*mu'))^2 ) for a dgCMatrix X
 * without additional storage
 */
// [[Rcpp::export(.icolSumSqdgC)]]
NumericVector colSumSqdgC(S4 X,NumericVector D,NumericVector mu)
{
  IntegerVector dims = X.slot("Dim");
  IntegerVector rowp = as<Rcpp::IntegerVector>(X.slot("i"));
  IntegerVector colp = as<Rcpp::IntegerVector>(X.slot("p"));
  NumericVector x = as<Rcpp::NumericVector>(X.slot("x"));

  int nrow = dims[0], ncol = dims[1];

  NumericVector v(ncol);
  double temp;
  for(int j=0;j<ncol;++j)
  {
    double vsum = 0.0;
    double muj = mu[j];
    int start = colp[j], end = colp[j+1];
    int iold = 0;
    for(int i=start; i<end;++i)
    {
      int irow = rowp[i];
      for(int l=iold;l<irow;++l)
      {
        temp = D[l]*muj;
        vsum += temp*temp;
      }
      temp = D[irow]*(x[i] - muj);
      vsum += temp*temp;
      iold = irow + 1;
    }
    for(int l=iold; l<nrow; ++l)
    {
      temp = D[l]*muj;
      vsum += temp*temp;
    }
    v[j] = 1.0/std::sqrt(vsum);
  }

  return(v);
}

/*
 * The following function computes the exact same thing as above but for regular matrices.
 *
 */
// [[Rcpp::export(.icolSumSq)]]
NumericVector colSumSq(NumericMatrix X,NumericVector D,NumericVector mu)
{
  int n = X.nrow(), p = X.ncol();
  NumericVector v(p);
  for(int j=0;j<p;++j)
  {
    double muj = mu[j];
    double vsum = 0.0;
    NumericMatrix::iterator ix = X.begin() + j*n;
    for(int i=0;i<n;++i)
    {
      double temp = D[i] * (*(ix++) - muj);
      vsum += temp*temp;
    }
    v[j] = 1/std::sqrt(vsum);
  }

  return(v);
}

// [[Rcpp::export(RXM)]]
NumericVector RXM_CC(NumericMatrix X,NumericVector ER)
{
  int n = X.nrow(), p = X.ncol();
  NumericVector v(p);
  
  for(int j=0;j<p;++j)
  {
    double sum = 0.0;
    NumericMatrix::iterator ix = X.begin() + j*n;
    for(int i=0;i<n;++i)
      sum += ER[i] * (*(ix++));
    
    v[j] = sum/n;
  }
  
  return(v);
}

//compute diagonal of cmat
// [[Rcpp::export(cmdg)]]
NumericVector cmdg_CC(NumericMatrix L, NumericVector D){
  int p = L.nrow(), q = L.ncol();
  NumericVector v(q);
  for(int k=0;k<q;++k){
    double  tmp = 1.0;
    for(int j=0;j<p;++j){
      tmp += 1/D[j]*L(j,k)*L(j,k);
    }
    v[k] = 1/tmp;
  }
  return(v);
}


// compute tao and iSxM
// [[Rcpp::export(taom)]]
List taom_CC(NumericMatrix X, NumericMatrix L, NumericVector D, NumericVector mu, NumericVector cmatdg)
{
  int n = X.nrow(), p = X.ncol(), q = L.ncol();
  NumericVector tao(n), xm(n), idm(p), summ2(q);
  NumericVector idx(p);
  
  for(int j=0;j<p;++j) {
	idm[j] = 1/D[j]*mu[j];
  }

  for(int k=0;k<q;++k){
    for(int j=0;j<p;++j){
      summ2[k] += idm[j]*L(j,k);
    }
  }
  
  for(int i=0;i<n;++i)
  {
    double sumt1 = 0.0, summ1 = 0.0, tempt = 0.0, tempm = 0.0;
    
    for(int k=0;k<q;++k){
      double sumt2 = 0.0;
      
      for(int j=0;j<p;++j){
        if(k == 0) {
          idx[j] = 1/D[j]*X(i,j);
          sumt1 += idx[j]*X(i,j);
          summ1 += idx[j]*mu[j];
          }
        sumt2 += idx[j]*L(j,k);
        }
        tempt += sumt2*sumt2*cmatdg[k];
        tempm += sumt2*summ2[k]*cmatdg[k];
      }
    tao[i] = sumt1 - tempt;
    xm[i] = summ1 - tempm;
  }
  
  return List::create(Named("tao") = tao,
                      Named("iSxM") = xm);
  
}

// compute only iSxM
// [[Rcpp::export(ism)]]
NumericVector ism_CC(NumericMatrix X, NumericMatrix L, NumericVector D, NumericVector mu, NumericVector cmatdg)
{
  int n = X.nrow(), p = X.ncol(), q = L.ncol();
  NumericVector xm(n), idm(p), summ2(q);
  NumericVector idx(p);
  
  for(int k=0;k<q;++k){
    for(int j=0;j<p;++j){
      if(k==0) {idm[j] = 1/D[j]*mu[j];}
      summ2[k] += idm[j]*L(j,k);
    }
  }
  
  for(int i=0;i<n;++i)
  {
    double summ1 = 0.0, tempm = 0.0;
    
    for(int k=0;k<q;++k){
      double sumt2 = 0.0;
      
      for(int j=0;j<p;++j){
        if(k == 0) {
          idx[j] = 1/D[j]*X(i,j);
          summ1 += idx[j]*mu[j];
        }
        sumt2 += idx[j]*L(j,k);
      }
      tempm += sumt2*summ2[k]*cmatdg[k];
    }
    xm[i] = summ1 - tempm;
  }
  
  return(xm);
  
}

// compute the SD vector correspond to the C
// [[Rcpp::export(mSD)]]
NumericVector mSD_CC(NumericMatrix X, NumericVector mu, NumericVector D,
                  NumericVector ER, NumericVector VR, NumericVector ybar )
    {
    int n = X.nrow(), p = X.ncol();
    NumericVector v(p);
    
    for(int j=0;j<p;++j)
    {
      double sum =0.0, temp = 0.0;
      for(int i=0;i<n;++i){
        sum += (ER[i]*ER[i]+VR[i])*X(i,j)*X(i,j)/n;
      }
      temp = sum - (2*ybar[j]*mu[j] - mu[j]*mu[j]);
      v[j] = std::sqrt(temp);
    }
    return(v);
    }

#endif

