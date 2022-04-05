#ifndef _SPH_CC_
#define _SPH_CC_
#define STRICT_R_HEADERS
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(sph)]]
void sph(SEXP X){
  
  SEXP dim     = Rf_getAttrib( X, R_DimSymbol );
  int nrow     = INTEGER(dim)[0];
  int ncol     = INTEGER(dim)[1];
  
  double * p = REAL(X) ;
  
  for(int i=0;i<nrow;++i)
  {
    double sum = 0.0;
    for(int j=0;j<ncol;++j){
      sum += p[i+nrow*j]*p[i+nrow*j]; //X(i,j) = p[i+nrow*j]
    }
    sum = std::sqrt(sum);
    
    for(int j=0;j<ncol;++j){
      p[i+nrow*j] /= sum;
    }
  }
}



#endif


