// disable all run-time checks in armadillo (e.g. bounds checking, etc.)
#define ARMA_NO_DEBUG

// we depend on the R package "RcppArmadillo"
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// we need OpenMP for parallelisation
// [[Rcpp::plugins(openmp)]]

// function for finding row medians (alternative to apply(depth, 1, median))
// requires integer type matrix as input, returns list of doubles
// [[Rcpp::export]]
std::vector<double> arma_rowMedians(const arma::imat &depth) {
    // number of rows
    const int nrows = depth.n_rows;

    // vector for storing the result
    std::vector<double> medians(nrows);

    // loop over the rows
    #pragma omp parallel for
    for (int i = 0; i < nrows; i++) {
        // convert the row to double type, to compute median correctly
        arma::rowvec row = arma::conv_to<arma::rowvec>::from(depth.row(i));

        // compute the median for this row
        medians[i] = arma::median(row);
    }

    return medians;
}

// function for finding row maximums (alternative to apply(mat, 1, max))
// requires integer type matrix as input, return list of integers
// [[Rcpp::export]]
std::vector<int> arma_rowMaximums(const arma::imat &mat) {
    const int nrows = mat.n_rows;

    // create vector to store the result
    std::vector<int> maximums(nrows);

    // loop over rows
    #pragma omp parallel for
    for (int i = 0; i < nrows; i++) {
        // find the maximum for this row
        maximums[i] = mat.row(i).max();
    }

    return maximums;
}

// C++ version of depth2K function
// [[Rcpp::export]]
Rcpp::NumericMatrix arma_depth2K(const Rcpp::NumericMatrix &A) {
    // create the output matrix (same size as input)
    Rcpp::NumericMatrix Aout(A.rows(), A.cols());

    // number of elements
    const long Asize = A.rows() * A.cols();

    // loop over elements in parallel and apply operation
    #pragma omp parallel for
    for (long i = 0; i < Asize; i++) {
        Aout[i] = 1.0 / pow(2.0, A[i]);
    }

    // return matrix
    return Aout;
}

// function for setting unused values of P0, P1 and genon01 to zero
// modifies the matrices in-place (i.e. doesn't return anything)
// [[Rcpp::export]]
void assignP0P1Genon01(Rcpp::NumericMatrix &P0, Rcpp::NumericMatrix &P1, Rcpp::NumericMatrix &genon01,
        const Rcpp::LogicalMatrix &usegeno, const Rcpp::NumericMatrix &dsub) {
    // number of elements (assumes all inputs are the same size!)
    const long size = P0.rows() * P0.cols();

    // loop over elements in parallel
    #pragma omp parallel for
    for (long i = 0; i < size; i++) {
        // set to zero if they match the conditions
        if (dsub[i] < 2.0) {
            P0[i] = 0.0;
            P1[i] = 0.0;
            genon01[i] = 0.0;
        }
        else if (!usegeno[i]) {
            P0[i] = 0.0;
            P1[i] = 0.0;
        }
    }

    // nothing to return, matrices are modified in-place
    return;
}
