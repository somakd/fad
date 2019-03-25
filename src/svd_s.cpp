
#ifndef _SVDS_XmD_CC
#define _SVDS_XmD_CC
#include <RcppEigen.h>
#include <Rcpp.h>
#include <SymEigs.h>
#include <GenEigs.h>
#include <RMatOp.h>
#include <R_ext/RS.h>
#include <Rinternals.h>

using namespace Spectra;


class XmD : public MatProd
{

private:
  const double *mat_ptr;
  const double  *mu_ptr;
  const double *D_ptr;
  double *vn;
  const int nrow;
  const int ncol;
  const double  BLAS_alpha;
  const int     BLAS_one;
public:

  XmD(SEXP mat_,SEXP mu_,SEXP D_, const int nrow_, const int ncol_) :
  mat_ptr(REAL(mat_)),
  mu_ptr(REAL(mu_)),
  D_ptr(REAL(D_)),
  nrow(nrow_),
  ncol(ncol_),
  BLAS_alpha(1.0),
  BLAS_one(1.0)
  {
    vn = new double[ncol];
  }

  int rows() const { return nrow; }
  int cols() const { return ncol; }

  // y_out = {(X - 1mu')*D} * x_in
  void perform_op(const double *x_in, double *y_out)
  {
    // compute Dx = D*x and set all y[i] = 1.0
    for(int i=0;i<ncol;++i) { vn[i] = D_ptr[i]*x_in[i];}
    for(int i=0;i<nrow;++i) y_out[i] = 1.0;

    double BLAS_beta = 0.0;
    for(int i=0;i<ncol;++i)
      BLAS_beta -= vn[i]*mu_ptr[i];

    // compute y = X*(Dx) + (-sum(mu*Dx))*y
    F77_CALL(dgemv)("N", &nrow, &ncol,
             &BLAS_alpha, mat_ptr, &nrow,
             (const double*)vn, &BLAS_one, &BLAS_beta,
             y_out, &BLAS_one);

  }

  // y = D*(X' - mu*1')*x
  void perform_tprod(const double *x_in, double *y_out)
  {
    double BLAS_beta = 0;
    for(int i=0;i<nrow;++i)
      BLAS_beta -= x_in[i];
    for(int i=0;i<ncol;++i) y_out[i] = mu_ptr[i];

    F77_CALL(dgemv)("T", &nrow, &ncol,
             &BLAS_alpha, mat_ptr, &nrow,
             x_in, &BLAS_one, &BLAS_beta,
             y_out, &BLAS_one);

    for(int i=0;i<ncol;++i)
      y_out[i] *= D_ptr[i];
  }


  ~XmD() {
    delete this->vn;
  }

};



// [[Rcpp::export(svds_XmD)]]
SEXP svds_XmD(SEXP X, SEXP mu,SEXP D, int k)
{
  BEGIN_RCPP
  SEXP dim     = Rf_getAttrib( X, R_DimSymbol );
  int m        = INTEGER(dim)[0];
  int n        = INTEGER(dim)[1];
  int nu       = k;
  int nv       = k;
  int ncv      = std::min((int)std::min(m,n),(int)std::max(2*k+1,21));
  double tol   = 1e-10;
  int maxitr   = 1000;


  MatProd *op_orig = new XmD(X,mu,D, m, n);
  MatProd *op;
  if(m > n)
    op = new SVDTallOp(op_orig);
  else
    op = new SVDWideOp(op_orig);

  SymEigsSolver<double, LARGEST_ALGE, MatProd> eigs(op, k, ncv);

  eigs.init();
  int nconv = eigs.compute(maxitr, tol);
  if(nconv < k)
    Rcpp::warning("only %d singular values converged, less than k = %d", nconv, k);
  nu = std::min(nu, nconv);
  nv = std::min(nv, nconv);

  Eigen::VectorXd evals = eigs.eigenvalues();
  Eigen::MatrixXd evecs = eigs.eigenvectors(std::max(nu, nv));

  Rcpp::NumericVector d(nconv);
  Rcpp::NumericMatrix u(m, nu), v(n, nv);
  int nops = 0;
  // Copy evals to d and take the square root
  std::copy(evals.data(), evals.data() + nconv, d.begin());
  for(int i=0;i<nconv;++i) d[i] = std::sqrt(d[i]);
  // Copy evecs to u or v according to the shape of A
  // If A is tall, copy evecs to v, otherwise copy to u
  if(m > n)
    std::copy(evecs.data(), evecs.data() + nv * n, v.begin());
  else
    std::copy(evecs.data(), evecs.data() + nu * m, u.begin());
  // Calculate the other one
  if(m > n)
  {
    // A = UDV', A'A = VD^2V', AV = UD, ui = A * vi / di
    // evecs has already been copied to v, so we can overwrite evecs
    for(int i = 0; i < nu; i++)
    {
      evecs.col(i) /= d[i];
      op_orig->perform_op(&evecs(0, i), &u(0, i));
      nops++;
    }
  } else {
    // A = UDV', AA' = UD^2U', A'U = VD, vi = A' * ui / di
    // evecs has already been copied to u, so we can overwrite evecs
    for(int i = 0; i < nv; i++)
    {
      evecs.col(i) /= d[i];
      op_orig->perform_tprod((const double *)&evecs(0, i), (double *)&v(0, i));
      nops++;
    }
  }
  Rcpp::RObject u_ret;
  Rcpp::RObject v_ret;

  if(nu > 0)  u_ret = u;  else  u_ret = R_NilValue;
  if(nv > 0)  v_ret = v;  else  v_ret = R_NilValue;

  delete op;
  delete op_orig;


  return Rcpp::List::create(
    Rcpp::Named("d")     = d,
    Rcpp::Named("u")     = u_ret,
    Rcpp::Named("v")     = v_ret,
    Rcpp::Named("niter") = eigs.num_iterations(),
    Rcpp::Named("nops")  = eigs.num_operations() * 2 + nops
  );
  END_RCPP
}


#endif

