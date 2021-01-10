
#ifndef _EIGS_SYM_RXmD_CC
#define _EIGS_SYM_RXmD_CC
#include <RcppEigen.h>
#include <Rcpp.h>
#include <SymEigs.h>
#include <GenEigs.h>
#include <RMatOp.h>
#include <R_ext/BLAS.h>  // for BLAS and F77_CALL

using namespace Spectra;

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

//input mat: X, vec: ER, VR, Ybar, mu, D
//--> DCD: C = ERX'ERX/N - muYbar' - Ybarmu' + mumu' + VRX'X/N
class RXmD : public MatProd
{

private:
  const double  *mat_ptr;
  const double  *mu_ptr;
  const double  *D_ptr;
  const double  *er_ptr;
  const double  *vr_ptr;
  const double  *ybar_ptr;
  double *v;
  const int nrow;
  const int ncol;
  const double BLAS_one;
  const int BLAS_one_int;
  const double BLAS_zero;
  const double BLAS_1byn;

public:

  RXmD(const NumericMatrix mat_,const  NumericVector mu_,const  NumericVector D_,
       const NumericVector er_, const NumericVector vr_, const NumericVector ybar_,
       const int nrow_, const int ncol_) :
  mat_ptr(REAL(mat_)),
  mu_ptr(REAL(mu_)),
  D_ptr(REAL(D_)),
  er_ptr(REAL(er_)),
  vr_ptr(REAL(vr_)),
  ybar_ptr(REAL(ybar_)),
  nrow(nrow_),
  ncol(ncol_),
  BLAS_one(1.0),
  BLAS_one_int(1),
  BLAS_zero(0.0),
  BLAS_1byn(1.0/nrow)
  {
    v = new double[nrow];
  }

  int rows() const { return ncol; }
  int cols() const { return ncol; }

  // compute DCD*x

  void perform_op(const double *x_in, double *y_out)
  {
    // compute y = D*x and set all v[i] = 0.0
    for(int i=0;i<ncol;++i) y_out[i] = D_ptr[i]*x_in[i];

    for(int i=0;i<nrow;++i) v[i] = 0;

    // compute v = X*y
    F77_CALL(dgemv)("N", &nrow, &ncol,
             &BLAS_one, mat_ptr, &nrow,
             (const double*)y_out, &BLAS_one_int, &BLAS_zero,
             v, &BLAS_one_int);

    double a = 0.0, b = 0.0;
    for(int j=0;j<ncol;++j)
    {
      a += mu_ptr[j]*y_out[j];
      b += ybar_ptr[j]*y_out[j];
    }
    b = a - b;

    // compute v = ER^2*v + V*v
    for(int i=0;i<nrow;++i)
      v[i] = (er_ptr[i]*er_ptr[i] + vr_ptr[i])*v[i];


    // compute y_out = b*mu - a*ybar
    for(int j=0;j<ncol;++j) y_out[j] = b*mu_ptr[j] - a*ybar_ptr[j];

    F77_CALL(dgemv)("T", &nrow, &ncol,
             &BLAS_1byn, mat_ptr, &nrow,
             (const double*)v, &BLAS_one_int, &BLAS_one,
             y_out, &BLAS_one_int);

    for(int i=0;i<ncol;++i)
      y_out[i] *= D_ptr[i];
    }

  void perform_tprod(const double *x_in, double *y_out) {
    perform_op(x_in,y_out);
  }



  ~RXmD() {
    delete this->v;
  }

};

// [[Rcpp::export]]
RcppExport SEXP eigs_sym_RXmD(
    SEXP Xmat, NumericVector mu,NumericVector D, NumericVector er, NumericVector vr,
    NumericVector ybar, int nev,int matclass)
{
  BEGIN_RCPP

  RXmD* op;
  int m,n;

  int maxitr = 1000;
  double tol = 1.0e-10;
  Rcpp::RObject evals, evecs;
  int nconv = 0, niter = 0, nops = 0;

  // check if X is matrix or not...
  if(matclass == 1) {
    SEXP dim = Rf_getAttrib( Xmat, R_DimSymbol );
    NumericMatrix X(Xmat);
    m = INTEGER(dim)[0];
    n = INTEGER(dim)[1];
    op = new RXmD(X,mu,D, er, vr, ybar, m, n);
  } else {
    Rcpp::stop("Currently only X of class matrix supported.\ndgCMatrix coming soon");
  }
  int ncv = std::min(std::min(n,m), std::max(2 * nev + 1, 20));

  SymEigsSolver<double, LARGEST_MAGN, RXmD> eigs(op, nev, ncv);

  eigs.init();
  nconv = eigs.compute(maxitr,tol);
  if(nconv < nev)
    Rcpp::warning("only %d eigenvalue(s) converged, less than k = %d",
                  nconv, nev);
  evals = Rcpp::wrap(eigs.eigenvalues());
  Rcpp::NumericMatrix evecs_ret(n, nconv);
  MapMat ma(evecs_ret.begin(), n, nconv);
  ma.noalias() = eigs.eigenvectors();
  evecs = evecs_ret;
  niter = eigs.num_iterations();
  nops = eigs.num_operations();

  delete op;

  return Rcpp::List::create(
    Rcpp::Named("values")  = evals,
    Rcpp::Named("vectors") = evecs,
    Rcpp::Named("nconv")   = nconv,
    Rcpp::Named("niter")   = niter,
    Rcpp::Named("nops")    = nops
  );

  END_RCPP
}

#endif
