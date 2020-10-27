// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <Rcpp.h>
#include <dqrng.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(dqrng)]]
// [[Rcpp::depends(ccaPP)]]

//' @export
//' @useDynLib chocR




//code copied from rcpp galery
//https://github.com/RcppCore/rcpp-gallery/blob/gh-pages/src/2013-07-13-dmvnorm_arma.Rmd
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
    arma::uword const n = trimat.n_cols;

    for(unsigned j = n; j-- > 0;){
        double tmp(0.);
        for(unsigned i = 0; i <= j; ++i)
            tmp += trimat.at(i, j) * x[i];
        x[j] = tmp;
    }
}

//code adapted from rcpp galery
//https://github.com/RcppCore/rcpp-gallery/blob/gh-pages/src/2013-07-13-dmvnorm_arma.Rmd
arma::vec dmvnrm_arma_fast(arma::mat const &x,
                           arma::rowvec const &mean,
                           arma::mat const &rooti,
                           double const rootisum,
                           double other_terms) {
    using arma::uword;
    uword const n = x.n_rows;
    arma::vec out(n);
    arma::rowvec z;
    for (uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);
        out(i) = other_terms - 0.5 * arma::dot(z, z);
    }
    return exp(out);
}



//code adapted from rcpp galery
//https://github.com/RcppCore/rcpp-gallery/blob/gh-pages/src/2013-07-13-dmvnorm_arma.Rmd
// [[Rcpp::export]]
arma::colvec dKernel(arma::mat const &grid,
                  arma::mat const &obs,
                  arma::vec const &probs,
                  arma::mat const &rooti) {
    double rootisum=arma::sum(log(rooti.diag()));
    double constants = -(double)grid.n_cols/2.0 * std::log(2.0 * M_PI);
    double other_terms = rootisum + constants;
    int nbx = grid.n_rows;
    int nnodes = obs.n_rows;
    arma::mat density(nbx, nnodes);
    for (int i=0; i<nnodes; ++i){
        density.col(i)=probs(i) * dmvnrm_arma_fast(grid, obs.row(i), rooti,
                    rootisum, other_terms);

    }
    return arma::sum(density, 1);
}

// [[Rcpp::export]]
arma::mat get_root_i(arma::mat const & chol){
    return arma::inv(trimatu(chol));
}

// [[Rcpp::export]]
arma::mat rKernel(int n, arma::mat x, arma::mat chol,
                  Rcpp::Nullable<Rcpp::NumericVector> probs = R_NilValue) {
    int ncols = chol.n_cols;
    arma::mat means(n, ncols);
    Rcpp::IntegerVector nodes=dqrng::dqsample_int(x.n_rows, n, true, probs);
    for (int i=0; i<n; ++i) {
        means.row(i) = x.row(nodes[i]);
    }

    arma::mat epsilon = arma::randn(n, ncols) * chol; //multivariate deviate

    return means + epsilon;
}

typedef void (*integrand) (unsigned ndim, const double *x, void *,
              unsigned fdim, double *fval);

double corKendall(arma:: vec const &vecx, arma::vec const &vecy, size_t n)
{
  arma::uvec indices = arma::sort_index(vecx);
  double * x = arma::vec(vecx.elem(indices)).memptr();
  double *y = arma::vec(vecy.elem(indices)).memptr();
  typedef double (*Fun)(double* arr1, double* arr2, size_t len, int cor) ;
    Fun fun = (Fun) R_GetCCallable( "pcaPP", "kendallNlogN" ) ;

    return fun(x, y, n, 1);
};

// [[Rcpp::export]]
arma::vec computeChoc(arma::mat const & grid,
                      Rcpp::List const &list_data,
                      Rcpp::List const &list_weights,
                      arma::mat const & rooti) {
    int ngrid=grid.n_rows;
    arma::vec result(ngrid);
    int nbyear=list_data.length();
    arma::vec years=arma::regspace(0,nbyear-1);
    arma::mat densities(ngrid, nbyear);
    for (int i=0; i<nbyear; ++i){
      arma::mat sdata = list_data[i];
      arma::vec weight = list_weights[i];
      densities.col(i) = dKernel(grid, sdata, weight, rooti);
    }
    for (int i=0; i< ngrid; ++i){
      arma::vec v (densities.n_cols);
      result(i) = corKendall(v,
             years, nbyear);
    }
    return result;
}



// [[Rcpp::export]]
arma::vec computeIvlev(arma::vec const & realised,
                       arma::vec const &avalaible) {
    return (realised - avalaible) / (realised + avalaible);
}





