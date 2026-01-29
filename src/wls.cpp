#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Implementation of the WLS loop
// Y: cells x detectors
// M: fluorophores x detectors
// background_noise: scalar
// MMt_inv: fluorophores x fluorophores (precomputed fallback)
// Returns A: cells x fluorophores

// [[Rcpp::export]]
arma::mat calc_wls_mat_cpp(arma::mat Y, arma::mat M, double background_noise, arma::mat MMt_inv) {
    int n_cells = Y.n_rows;
    int n_fluor = M.n_rows;
    int n_detectors = M.n_cols;

    // Check dimensions
    if (Y.n_cols != n_detectors) {
        stop("Y and M must have same number of columns (detectors)");
    }

    arma::mat A(n_cells, n_fluor, fill::zeros);
    arma::mat Mt = M.t();
    arma::rowvec zeros_det = arma::zeros<arma::rowvec>(n_detectors);

    for (int i = 0; i < n_cells; ++i) {
        arma::rowvec yi = Y.row(i);

        // Calculate weights
        // R: weights_i <- 1 / (pmax(Y[i, ], 0) + background_noise)
        arma::rowvec weights = 1.0 / (arma::max(yi, zeros_det) + background_noise);

        // Construct weighted matrix
        // R: M %*% Wi %*% Mt
        // Wi is diag(weights).
        // In Armadillo: (M.each_row() % weights) * Mt
        // Note: weights is 1 x D. M is F x D.
        // We want to scale each column j of M by weights[j].
        // M.each_row() % weights multiplies each row of M by weights.
        // M row k element j is M_kj. Result is M_kj * weights_j.
        // This effectively scales column j by weights_j.

        arma::mat Mw = M;
        Mw.each_row() %= weights;

        arma::mat MWMt_i = Mw * Mt;

        // Check condition
        // R: if (rcond(MWMt_i) > 1e-10)
        double rcond_val = arma::rcond(MWMt_i);

        if (rcond_val > 1e-10) {
            // A[i, ] <- Y[i, ] %*% Wi %*% Mt %*% solve(MWMt_i)
            // Y[i, ] %*% Wi is yi % weights (element-wise)
            arma::rowvec yw = yi % weights;

            // yw %*% Mt
            arma::rowvec rhs = yw * Mt;

            // rhs %*% solve(MWMt_i)
            // We need to solve x * MWMt_i = rhs -> MWMt_i^T * x^T = rhs^T -> MWMt_i * x^T = rhs^T (since symm)
            // x^T = solve(MWMt_i, rhs^T)

            arma::vec sol = arma::solve(MWMt_i, rhs.t());
            A.row(i) = sol.t();
        } else {
            // A[i, ] <- Y[i, ] %*% Mt %*% MMt_inv
            // Fallback to OLS
            A.row(i) = yi * Mt * MMt_inv;
        }
    }

    return A;
}

extern "C" SEXP _spectrQC_calc_wls_mat_cpp(SEXP YSEXP, SEXP MSEXP, SEXP bgSEXP, SEXP MMt_invSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type background_noise(bgSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type MMt_inv(MMt_invSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_wls_mat_cpp(Y, M, background_noise, MMt_inv));
    return rcpp_result_gen;
    END_RCPP
}
