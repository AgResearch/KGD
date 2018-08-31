// disable all run-time checks in armadillo (e.g. bounds checking, etc.)
#define ARMA_NO_DEBUG

// we depend on the R package "RcppArmadillo"
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// we need OpenMP for parallelisation
#include <omp.h>
// [[Rcpp::plugins(openmp)]]


// helper function for deciding number of threads
static int check_nThreads(int nThreads) {
	// maximum number of threads available
	int maxThreads = omp_get_max_threads();

	if (nThreads <= 0) {
		// if nThreads is set to zero then use everything
		nThreads = maxThreads;
	}
	else if (nThreads > maxThreads) {
		// don't allow more threads than the maximum available
		nThreads = maxThreads;
	}

	return nThreads;
}


// function for finding row medians (alternative to apply(depth, 1, median))
// requires integer type matrix as input, returns list of doubles
// [[Rcpp::export]]
std::vector<double> rcpp_rowMedians(const arma::imat &depth, int nThreads) {
	// set up number of threads
	nThreads = check_nThreads(nThreads);

    // number of rows
    const int nrows = depth.n_rows;

    // vector for storing the result
    std::vector<double> medians(nrows);

    // loop over the rows
    #pragma omp parallel for num_threads(nThreads)
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
std::vector<int> rcpp_rowMaximums(const arma::imat &mat, int nThreads) {
	// set up number of threads
	nThreads = check_nThreads(nThreads);

	// number of rows
    const int nrows = mat.n_rows;

    // create vector to store the result
    std::vector<int> maximums(nrows);

    // loop over rows
    #pragma omp parallel for num_threads(nThreads)
    for (int i = 0; i < nrows; i++) {
        // find the maximum for this row
        maximums[i] = mat.row(i).max();
    }

    return maximums;
}

// C++ version of depth2K function
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_depth2K(const Rcpp::NumericMatrix &A, int nThreads) {
	// set up number of threads
	nThreads = check_nThreads(nThreads);

    // create the output matrix (same size as input)
    Rcpp::NumericMatrix Aout(A.rows(), A.cols());

    // number of elements
    const long Asize = A.rows() * A.cols();

    // loop over elements in parallel and apply operation
    #pragma omp parallel for num_threads(nThreads)
    for (long i = 0; i < Asize; i++) {
        Aout[i] = 1.0 / pow(2.0, A[i]);
    }

    // return matrix
    return Aout;
}

// C++ version of depth2Kmodp function
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_depth2Kmodp(const Rcpp::NumericMatrix &depthvals, double modp, int nThreads) {
	// set up number of threads
	nThreads = check_nThreads(nThreads);

    // create matrix for storing the result
    Rcpp::NumericMatrix result(depthvals.rows(), depthvals.cols());

    // size of the matrix
    const long size = depthvals.rows() * depthvals.cols();

    // loop over the elements in parallel
    #pragma omp parallel for num_threads(nThreads)
    for (long i = 0; i < size; i++) {
        double value = 0.5 * pow(modp, depthvals[i] - 1.0);
        result[i] = (value == 0) ? 1.0 : value;
    }
    return result;
}

// C++ version of depth2Kbb function
// [[Rcpp::export]]
    Rcpp::NumericMatrix rcpp_depth2Kbb(const Rcpp::NumericMatrix & depthvals, int nThreads, const double alph = 9999) {
        // set up number of threads
        nThreads = check_nThreads(nThreads);
        // create matrix for storing the result
        Rcpp::NumericMatrix result(depthvals.rows(), depthvals.cols());
        // size of the matrix
        const long size = depthvals.rows() * depthvals.cols();
        // precompute factor
        const double factor = 1.0/R::beta(alph, alph);
        // loop over the elements in parallel
        #pragma omp parallel for num_threads(nThreads)
        for (long i = 0; i < size; i++) {
            result[i] = R::beta(alph, depthvals[i] + alph) * factor;
        }
        return result;
    }


// function for setting unused values of P0, P1 and genon01 to zero
// modifies the matrices in-place (i.e. doesn't return anything)
// [[Rcpp::export]]
void rcpp_assignP0P1Genon01(Rcpp::NumericMatrix &P0, Rcpp::NumericMatrix &P1, Rcpp::NumericMatrix &genon01,
        const Rcpp::LogicalMatrix &usegeno, const Rcpp::NumericMatrix &dsub, int nThreads) {
	// set up number of threads
	nThreads = check_nThreads(nThreads);

    // number of elements (assumes all inputs are the same size!)
    const long size = P0.rows() * P0.cols();

    // loop over elements in parallel
    #pragma omp parallel for num_threads(nThreads)
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
