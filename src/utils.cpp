#include <RcppArmadillo.h>
#include "brent.hpp"
#include <functional>


using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

////////////////////////////////////////
// Auxiliary functions
////////////////////////////////////////


arma::vec removeElement(vec x, int ix){
  // Auxiliary function.
  // Given a vector x, returns a vector out without the element with index ix.
  // Note: this is 0-based indexing like in C++ and not 1-based indexing like in R
  int n = x.n_elem;
  vec out;
  // If ix is not the first or last element
  if(ix != 0 & ix != n-1){
    out = join_vert(x.subvec(0, ix-1), x.subvec(ix + 1, n-1));
  } else{
    if(ix == 0){
      // If ix is the first element
      out = x.subvec(1, n-1);
    } else{
      // If ix is the last element
      out = x.subvec(0, n-2);
    }
  }
  return(out);
}



// [[Rcpp::export]]
arma::mat computeCoxDirection(arma::vec& x, int comp, int n_points, int verbose){
  // Function that returns a discretization of the Cox direction for vector x in component comp.
  // Returns a matrix of dimension (cox_dir_n_elems, q), where cox_dir_n_elems is roughly equal to n_points and q is the length of x.
  // Input:
  //     x: q dimensional vector of proportions. Must sum up to 1.
  //        There's no input check of this because this operation is done many times in the coordinate exchange algorithm and
  //         because the x vector is fed by another function.
  //     comp: component in which the Cox direction is computed. Must be an integer between 1 and q.
  //     n_points: Number of points to use in the discretization.
  //     verbose: integer that expresses the level of verbosity. Mainly used in other functions and too much useful by itself.

  if(verbose >= 2) Rcout << "Computing Cox direction" << std::endl;

  int i = comp - 1;
  int q = x.n_elem;

  //   points in Cox's direction
  vec seq_points = linspace<vec>(0, 1, n_points);

  vec diffs = seq_points - x(i);

  int ix1 = index_min(abs(diffs));

  // Auxiliary vector 1
  vec cox_direction_aux_1 = seq_points - diffs(ix1);

  // Auxiliary vector 2
  // copy all elements of cox_direction_aux_1 except the one with index ix1
  // the equivalent in R: cox_direction_aux_2 = cox_direction_aux_1[-ix1]
  vec cox_direction_aux_2 = removeElement(cox_direction_aux_1, ix1);

  uvec bigger_zero = find(cox_direction_aux_2 >= 0); // indices of elements bigger than 0
  uvec less_one = find(cox_direction_aux_2 <= 1); // indices of elements less than 1
  uvec betw_0_1 = intersect(bigger_zero, less_one); // indices of elements between 0 and 1

  // Auxiliary vector 3
  vec cox_direction_aux_3 = join_vert(zeros(1), cox_direction_aux_2.elem(betw_0_1), ones(1));

  // Number of elements in final vector
  int cox_dir_n_elems = cox_direction_aux_3.n_elem;

  // Initialize final matrix with random numbers
  arma::mat cox_direction = arma::mat(cox_dir_n_elems, q, fill::randu);

  // Fill corresponding column
  cox_direction.col(i) = cox_direction_aux_3;

  // Vector of differences
  vec deltas = cox_direction_aux_3 - x(i);

  // Fill matrix
  for(int n = 0; n < cox_dir_n_elems; n++){

    // recompute proportions:
    vec setDiff_aux = linspace<vec>(0, q-1, q);
    vec setDiff = removeElement(setDiff_aux, i);
    int j;
    double result;

    // Iterate over ingredient proportions
    for(int j_aux = 0; j_aux < setDiff.n_elem; j_aux++){
      j = setDiff(j_aux);
      if(abs(1 - x(i)) < 1e-16) {
        // In case x(i) is numerically 1
        result = (1 - cox_direction(n, i))/(q-1);
      } else{
        // In case x(i) is not numerically 1
        result = x(j) - deltas(n)*x(j)/(1 - x(i));
      }
      cox_direction(n, j) = result;
      j++;
    } // end for j

    // Check that the computed directions are not smaller than 0 or bigger than 1 because of numerical error
    if(any(cox_direction.row(n) < -1e-10 || cox_direction.row(n) > 1 + 1e10)) {
      Rcout << "Value out of bounds while computing Cox direction." << std::endl;
      Rcout << "The error ocurred with the following vector in component " << comp << std::endl;
      Rcout << x << std::endl;
      Rcout << "Cox direction computed as:" << std::endl;
      Rcout << cox_direction << std::endl;


      for(int j = 0; j < q; j++){
        cox_direction(n, j) = 1.0/q;
      }

      Rcout << "Filling Cox direction with 1/q. Corrected Cox direction:" << std::endl;
      Rcout << cox_direction << std::endl << std::endl << std::endl;

      warning("Value out of bounds while computing Cox direction. Filling with 1/q.\n");

    }
  } // end for n
  return(cox_direction);
}





////////////////////////////////////////
// Gaussian model
////////////////////////////////////////



arma::mat getScheffeGaussianOrder2(arma::mat& X){
  int q = X.n_cols;
  int n = X.n_rows;
  int n_col_X_m = q + (q-1)*q/2;
  arma::mat X_m(n, n_col_X_m);

  // Copy X matrix into first q columns of X_m
  for(int i = 0; i < n; i++){
    for(int j = 0; j < q; j++){
      X_m(i,j) = X(i,j);
    }
  }

  // Fill rest of matrix element-wise
  int k = q-1;
  for(int i = 0; i < (q-1); i++){
    for(int j = i+1; j < q; j++){
      k = k+1;
      for(int row = 0; row < n; row++){
        X_m(row,k) = X(row,i)*X(row,j);
      }
    }
  }

  return X_m;
}



arma::mat getScheffeGaussianOrder3(arma::mat& X){
  int q = X.n_cols;
  int n = X.n_rows;

  // This probably slows down the whole thing
  arma::mat X_ord_2 = getScheffeGaussianOrder2(X);

  int n_col_X_ord_2 = X_ord_2.n_cols;

  int n_col_X_m = X_ord_2.n_cols;
  // compute number of columns in X_m
  // There's a formula to find this number, but I'll work it out later.
  for(int i = 0; i < (q-2); i++){
    for(int j = (i+1); j < (q-1); j++){
      for(int k = (j+1); k < q; k++){
        n_col_X_m++;
      }
    }
  }

  arma::mat X_m(n, n_col_X_m);

  // Copy X_ord_2 matrix into first q columns of X_m
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n_col_X_ord_2; j++){
      X_m(i,j) = X_ord_2(i,j);
    }
  }

  // // Fill rest of matrix column-wise
  int l = q - 1;
  for(int i = 0; i < (q-1); i++){
    for(int j = i+1; j < q; j++){
      l++;
      for(int row = 0; row < n; row++){
        X_m(row,l) = X(row,i)*X(row,j);
      }
    }
  }

  l = X_ord_2.n_cols - 1;
  for(int i = 0; i < (q - 2); i++){
    for(int j = i + 1; j < (q - 1); j++){
      for(int k = j + 1; k < q; k++){
        l++;
        for(int row = 0; row < n; row++){
          X_m(row,l) = X(row,i)*X(row,j)*X(row,k);
        }
      }
    }
  }


  return X_m;
}




// [[Rcpp::export]]
arma::mat getScheffeGaussian(arma::mat& X, int order){

  // not the most elegant, but works
  bool flag = (order != 1 & order != 2 & order != 3);

  if(flag){
    stop("Inadmissible value for order. Must be 1, 2 or 3");
  }

  arma::mat X_m;

  if(order == 1)  X_m = X;
  else{
    if(order == 2){
      X_m = getScheffeGaussianOrder2(X);
    } else {
      X_m = getScheffeGaussianOrder3(X);
    }
  }
  // arma::mat_out = as<mat>(X_m);
  // return as<mat>(X_m);
  return X_m;
}



// [[Rcpp::export]]
double getDCritValueGaussian(arma::mat& X, int order){
  arma::mat X_m = getScheffeGaussian(X, order);
  arma::mat X_mT = trans(X_m);
  arma::mat I = X_mT * X_m; // Information matrix
  double eff_crit; // We want to minimize this

  // Attempt to do a Cholesky decomposition on the information matrix
  arma::mat L;
  double log_det_I;
  try{
    L = chol(I);
    // Compute the determinant of information matrix using the decomposition
    log_det_I = 2*sum(log(L.diag()));
    eff_crit = -log_det_I/X_m.n_cols + log(X_m.n_rows);
  }
  catch(const std::runtime_error& e){
    // If Cholesky decomposition fails, it is likely because information matrix
    // was not numerically positive definite.
    // If this happens, it is probably because a numerical inestability.
    // The function then returns the efficiency value as a big positive number, this
    // way the algorithm does nothing in this iteration because the algorithm thinks
    // there was no improvement when swapping the proportions.
    eff_crit = 100;
    Rcpp::warning("Error in Cholesky decomposition with message: ", e.what(), "\nReturning eff_crit = ", eff_crit);
  }

  return eff_crit;
}




// [[Rcpp::export]]
double getICritValueGaussian(arma::mat& X, int order, int q, arma::mat& W){

  if(order != 3) stop("Only special cubic models are allowed (i.e., order = 3).");

  arma::mat X_m = getScheffeGaussian(X, order);
  arma::mat X_mT = trans(X_m);
  arma::mat I = X_mT * X_m; // Information matrix

  double eff_crit; // We want to minimize this


  // Attempt to do a Cholesky decomposition on the information matrix
  arma::mat L;
  try{
    L = chol(I);

    arma::mat A = solve(trimatl(L.t()), W);
    arma::mat C = solve(trimatu(L), A);
    eff_crit = trace(C);
  }
  catch(const std::runtime_error& e){
    // If Cholesky decomposition fails, it is likely because information matrix
    // was not numerically positive definite.
    // If this happens, it is probably because a numerical inestability.
    // The function then returns the efficiency value as a big positive number, this
    // way the algorithm does nothing in this iteration because the algorithm thinks
    // there was no improvement when swapping the proportions.
    eff_crit = 100;
    Rcpp::warning("Error in Cholesky decomposition with message: ", e.what(), "\nReturning eff_crit = ", eff_crit);
  }

  return eff_crit;
}





// [[Rcpp::export]]
double getOptCritValueGaussian(arma::mat& X, int order, int q, int opt_crit, arma::mat& W){
  //  opt_crit: optimality criterion: 0 (D-optimality) or 1 (I-optimality)
  if(opt_crit == 0){
    return(getDCritValueGaussian(X, order));
  } else{
    return(getICritValueGaussian(X, order, q, W));
  }
}




arma::mat findBestCoxDirGaussianDiscrete(
    arma::mat& cox_dir, arma::mat& X_in, int k, int order, double opt_crit_value_best,
    int opt_crit, arma::mat& W) {
  arma::mat X = X_in;
  int n_col_X = X.n_cols;
  arma::vec x_k(n_col_X);

  double opt_crit_value_j;
  int n_cox_points = cox_dir.n_rows;
  for(int j = 0; j < n_cox_points; j++){

    // In Rcpp: x_k = X(k-1,_);
    for(int elem = 0; elem < n_col_X; elem++){
      x_k(elem) = X(k-1, elem);
    }

    // In Rcpp: X(k-1,_) = cox_dir(j, _);
    for(int elem = 0; elem < n_col_X; elem++){
      X(k-1, elem) = cox_dir(j, elem);
    }

    opt_crit_value_j = getOptCritValueGaussian(X, order, n_col_X, opt_crit, W);

    // If new optimality criterion value is better, then keep the new one.
    if(opt_crit_value_j < opt_crit_value_best) {
      // This design has a better optimality criterion value, so we keep the design and update the best value
      opt_crit_value_best = opt_crit_value_j;
    } else{
      // This design does not have a better optimality criterion value, so we return to the old design.
      // In Rcpp: X(k-1,_) = x_k;
      for(int elem = 0; elem < n_col_X; elem++){
        X(k-1,elem) = x_k(elem);
      }

    }
  }

  return X;
}






// [[Rcpp::export]]
arma::mat changeIngredientDesignGaussian(double theta, arma::mat& X, int i, int j){
  // Returns a new design matrix Y changing the j-th ingredient in i-th observation is changed to theta
  // in the original design matrix X.
  // theta must be between 0 and 1 because it's an ingredient proportion.
  // j and i are 0-indexed.


  // Create new matrix Y that is identical to the one pointed by X.
  // Note: This is the easiest way to do it because we have to modify a row in this matrix.
  // A more computationally effective way would be to only store the new modified vector since
  // we don't need a copy of the whole matrix. But to do that I would have to either modify some
  // existing functions, or create some new ones, or both. IDK if the gain in performance is worth it.
  arma::mat Y = X;
  int q = Y.n_cols;

  // Create a vector with the i-th row of the design matrix X
  arma::vec x_i(q);
  for(int col = 0; col < q; col++){
    x_i(col) = X(i, col);
  }

  // delta = theta - x_i[i]
  double delta = theta - x_i(j);

  // recompute proportions:
  vec setDiff_aux = linspace<vec>(0, q-1, q);
  vec setDiff = removeElement(setDiff_aux, j);
  int k;
  double result;

  for(int k_aux = 0; k_aux < setDiff.n_elem; k_aux++){
    k = setDiff(k_aux);

    if(abs(1 - x_i(j)) < 1e-16) {
      // In case x_i(j) is numerically 1, it will return a numeric zero such that the vector sums up to 1
      // Almost the same as doing result = 0;
      result = (1 - x_i(j))/(q-1);
    } else{
      // Other case
      result = x_i(k) - delta*x_i(k)/(1 - x_i(j));
    }

    x_i(k) = result;
  }

  x_i(j) = theta;

  // replace x_i row with the recomputed proportions according to Cox direction
  for(int col = 0; col < q; col++){
    Y(i, col) = x_i(col);
  }

  return(Y);

}




// [[Rcpp::export]]
double efficiencyCoxScheffeGaussian(double theta, arma::mat& X, int i, int j, int order,
                                    int opt_crit, arma::mat& W){
  // Computes efficiency criterion of a design matrix X but where the j-th ingredient in i-th observation is changed to theta.
  // theta must be between 0 and 1 because it's an ingredient proportion.
  // j and i are 0-indexed.
  // We want to minimize this.

  arma::mat Y = changeIngredientDesignGaussian(theta, X, i, j);
  int q = Y.n_cols;

  // Return utility function value. We want to minimize this.
  return(getOptCritValueGaussian(Y, order, q, opt_crit, W));
}



// Help from:
// https://stackoverflow.com/questions/34994475/rcpp-function-to-construct-a-function
// https://codereview.stackexchange.com/questions/103762/implementation-of-brents-algorithm-to-find-roots-of-a-polynomial

// Thanks to Norma and Manuel who helped me out with lambda functions in C++:
// https://github.com/NormaVTH
// https://github.com/manuelalcantara52




arma::mat findBestCoxDirGaussianBrent(
    arma::mat& X, int i, int j, int order, int opt_crit, arma::mat& W,
    double lower = 0, double upper = 1, double tol = 0.0001) {

  auto f = [&X, i, j, order, opt_crit, &W](double theta){
    return efficiencyCoxScheffeGaussian(theta, X, i, j, order, opt_crit, W);
  };

  double theta_brent, theta_star, f_star;
  double f_brent = brent::local_min_mb(lower, upper, tol, f, theta_brent);
  // Check end points
  double f_lower = efficiencyCoxScheffeGaussian(lower, X, i, j, order, opt_crit, W);
  double f_upper = efficiencyCoxScheffeGaussian(upper, X, i, j, order, opt_crit, W);

  f_star = std::min({f_lower, f_upper, f_brent});

  if(f_star == f_brent){
    theta_star = theta_brent;
  } else{
    if(f_star == f_lower){
      theta_star = lower;
    } else{
      theta_star = upper;
    }
  }

  return(changeIngredientDesignGaussian(theta_star, X, i, j));

}






// [[Rcpp::export]]
Rcpp::List mixtureCoordinateExchangeGaussian(
    arma::mat X_orig, int order, int max_it, int verbose, int opt_crit,
    arma::mat W, int opt_method, double lower, double upper, double tol, int n_cox_points){
  // Performs the coordinate exchange algorithm for a Multinomial Logit Scheffé model.
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
  // X: armadillo matrix with dimensions (n, q) where:
  //    q is the number of ingredient proportions
  //    n is the number of runs
  // n_cox_points: Number of points to use in the discretization of Cox direction.
  // max_it: Maximum number of iteration that the coordinate exchange algorithm will do.
  // verbose: level of verbosity. 6 levels, in which level prints the previous plus additional things:
  //    1: Print the optimality criterion value in each iteration and a final summary
  //    2: Print the values of k, s, i, and optimality criterion value in each subiteration
  //    3: Print the resulting X after each iteration, i.e., after each complete pass on the data
  //    4: Print optimality criterion value for each point in the Cox direction discretization
  //    5: Print the resulting X and information matrix after each subiteration
  //    6: Print the resulting X or each point in the Cox direction discretization
  // opt_crit: optimality criterion: 0 (D-optimality) or 1 (I-optimality)
  // opt_method: optimization method (0 for Brent, 1 for discrete approximation).
  // W: moment matrix for I-optimality
  //
  // Returns an Rcpp::List object with the following objects:
  //    X_orig: The original design. Armadillo matrix with dimensions (n, q).
  //    X: The optimized design. Armadillo matrix with dimensions (n, q).
  //    opt_crit_value_orig: optimality criterion value of the original design.
  //    opt_crit_value: optimality criterion value of the optimized design.
  //    n_iter: Number of iterations performed.

  // Does not do input checks because the R wrapper function does them.

  // Create new matrix, otherwise it is modified in R too
  arma::mat X = X_orig;

  int n_runs = X.n_rows;
  int q = X.n_cols;

  // Create matrix with appropriate dimensions for Cox direction in each iteration
  arma::mat cox_dir(n_cox_points, q);

  // Vector of ingredient proportions
  arma::vec x(q);

  double opt_crit_value_orig = getOptCritValueGaussian(X, order, q, opt_crit, W);
  double opt_crit_value_best = opt_crit_value_orig;
  double opt_crit_value_aux = -1e308; // -Inf

  // Coordinate exchanges
  int it = 0;
  while(it < max_it){
    it = it + 1;
    if(verbose >= 1) Rcout << "Iter: " << it << ", Optimality criterion value: " << opt_crit_value_best << std::endl;

    // If there was no improvement in this iteration
    if(abs(opt_crit_value_aux - opt_crit_value_best) < 1e-16) break;

    opt_crit_value_aux = opt_crit_value_best;

    for(int k = 1; k <= n_runs; k++){

      // Checking interruption every 2 iterations
      // As in https://teuder.github.io/rcpp4everyone_en/270_miscellaneous.html
      if (k % 2 == 0){
        Rcpp::checkUserInterrupt();
      }

      for(int i = 0; i < q; i++){

        if(verbose >= 2) Rcout << "\nIter: " << it <<  ", k = " << k << ", i = " << i << std::endl;

        // populate x vector with corresponding ingredient proportions
        for(int l = 0; l < q; l++){
          x(l) = X(k-1, l);
        }

        if(opt_method == 0){
          // arma::mat& X, int i, int j, int order, int opt_crit, arma::mat& W,
          // double lower = 0, double upper = 1, double tol = 0.0001
          X = findBestCoxDirGaussianBrent(X, k-1, i, order, opt_crit, W, lower, upper, tol);

        } else{
          cox_dir = computeCoxDirection(x, i+1, n_cox_points, verbose);
          X = findBestCoxDirGaussianDiscrete(cox_dir, X, k, order, opt_crit_value_best, opt_crit, W);
        }


        opt_crit_value_best = getOptCritValueGaussian(X, order, q, opt_crit, W);

        if(verbose >= 2) Rcout << "Opt-crit-value: " << opt_crit_value_best << std::endl;

        if(verbose >= 5){
          Rcout << "X =\n" << X << std::endl;
        }

      } // end for i

    } // end for k

    if(verbose >= 3) Rcout << "X =\n" << X << std::endl;

    if(verbose >= 2) Rcout << std::endl << std::endl;

  } // end while

  if(verbose >= 1){
    Rcout << std::endl;
    Rcout << "Original Optimality criterion value: " << opt_crit_value_orig;
    Rcout << std::endl;
    Rcout << "Final Optimality criterion value: " << opt_crit_value_best;
    Rcout << std::endl;
    Rcout << "Number of iterations: " << it;
    Rcout << std::endl;
  }

  // return object
  return Rcpp::List::create(
    _["X_orig"] = X_orig,
    _["X"] = X,
    _["opt_crit_value_orig"] = opt_crit_value_orig,
    _["opt_crit_value"] = opt_crit_value_best,
    _["n_iter"] = it,
    _["opt_crit"] = opt_crit
  );

} // end function

















////////////////////////////////////////
// Multinomial logit (MNL) model
////////////////////////////////////////



// [[Rcpp::export]]
arma::mat getXsMNL(arma::cube& X, int s){
  // Function that returns the design matrix of choice set s.
  // Final matrix is of dimension (J, m-1) with m = (q^3 + 5*q)/6
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
  // Input should be a design cube X of dimensions (q, J, S) and integer s, corresponding to a choice set in 1 to S.
  int q = X.n_rows;
  int J = X.n_cols;
  int m = (q*q*q + 5*q)/6;

  // Initialize array with zeros
  arma::mat Xs(J, m-1, fill::zeros);

  // Column counter. Equation has three terms, so the final matrix is populated in three parts.
  int col_counter = 0;

  // First part
  for(int i = 0; i < q-1; i++){
    for(int j = 0; j < J; j++){
      // subtract 1 from s because it is 1-indexed
      Xs(j, col_counter) = X(i, j, s-1);
    }
    col_counter++;
  }

  // second part
  for(int i = 0; i < q-1; i++){
    for(int k = i+1; k < q; k++){
      for(int j = 0; j < J; j++){
        Xs(j, col_counter) = X(i, j, s-1)*X(k, j, s-1);
      }
      col_counter++;
    }
  }

  // third part
  for(int i = 0; i < q-2; i++){
    for(int k = i+1; k < q-1; k++){
      for(int l = k+1; l < q; l++){
        for(int j = 0; j < J; j++){
          Xs(j, col_counter) = X(i, j, s-1)*X(k, j, s-1)*X(l, j, s-1);
        }
        col_counter++;
      }
    }
  }

  return Xs;
}




// [[Rcpp::export]]
arma::vec getUsMNL(arma::cube& X, arma::vec& beta, int s, arma::mat& Xs){
  // Function that returns the utility vector of choice set s.
  // Final vector is of length J.
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
  // Input:
  //     X: design cube of dimensions (q, J, S)
  //     beta: parameter vector. Must be of length m, with m = (q^3 + 5*q)/6
  //     s: integer s, corresponding to a choice set in 1 to S.
  //     Xs: design matrix of choice set s. Must be of dimension (J, m-1), with m = (q^3 + 5*q)/6

  int J = X.n_cols;
  int q = X.n_rows;

  int m = Xs.n_cols + 1; // m = (q^3 + 5*q)/6

  // Check input dimensions
  if(m != beta.n_elem) stop("Incompatible q in beta and X");

  // Create auxiliary vector
  arma::vec beta2(m-1);

  // compute beta_i_star = beta_i - beta_q
  for(int i = 0; i < q-1; i++){
    beta2(i) = beta(i) - beta(q-1);
  }

  for(int i = q-1; i < m-1; i++){
    beta2(i) = beta(i+1);
  }

  arma::vec Us(J);

  Us = Xs*beta2;
  return Us;
}


// [[Rcpp::export]]
arma::vec getPsMNL(arma::cube& X, arma::vec& beta, int s, arma::mat& Xs){
  // Function that returns the probability vector of choice set s, based on the softmax function.
  // Final vector is of length J.
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
  // Input:
  //     X: design cube of dimensions (q, J, S)
  //     beta: parameter vector. Must be of length m, with m = (q^3 + 5*q)/6
  //     s: integer s, corresponding to a choice set in 1 to S.
  //     Xs: design matrix of choice set s. Must be of dimension (J, m-1), with m = (q^3 + 5*q)/6

  int J = X.n_cols;

  arma::vec Us(J);
  Us = getUsMNL(X, beta, s, Xs);

  arma::vec exp_Ujs(J);
  arma::vec P(J);

  // subtracting the maximum value to avoid numerical overflow
  exp_Ujs = exp(Us - max(Us));

  double sum_exp_Ujs = sum(exp_Ujs);

  P = exp_Ujs/sum_exp_Ujs;

  // Check for numerical inestabilities
  if(abs(sum(P) - 1) > 1e-10) warning("Sum may not be numerically equal to 1.");

  return P;
}





// [[Rcpp::export]]
arma::mat getInformationMatrixMNL(arma::cube& X, arma::vec& beta){
  // Function that returns the information matrix for design cube X and parameter vector beta.
  // It is the sum of the information matrices of the S choice sets.
  // Final matrix is of dimension (m-1, m-1), with m = (q^3 + 5*q)/6
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
  // Input:
  //     X: design cube of dimensions (q, J, S)
  //     beta: parameter vector. Must be of length m, with m = (q^3 + 5*q)/6

  int J = X.n_cols;
  int S = X.n_elem/(X.n_cols*X.n_rows);
  int m = beta.n_elem;

  arma::mat Xs(J, m-1);
  arma::mat I(m-1, m-1, fill::zeros);
  arma::mat identity(J, J, fill::eye);
  arma::mat ps_ps_t(J, J);
  arma::mat middle(J, J);
  arma::vec ps;

  // Compute information matrix for each choice set s, and sum.
  for(int s = 1; s <= S; s++){
    Xs = getXsMNL(X, s);
    ps = getPsMNL(X, beta, s, Xs);;
    ps_ps_t = ps*ps.t();
    middle = ps_ps_t;
    middle.diag() = ps_ps_t.diag() - ps;
    I = I - (Xs.t())*middle*Xs;
  }
  return I;
}







// // [[Rcpp::export]]
// double getDCritValueMNL(arma::cube& X, arma::mat& beta_mat, int verbose){
//   // Function that returns the D criterion value for design cube X and parameter vector beta.
//   // The D-optimality criterion seeks to maximize the determinant of the information matrix or
//   //   minimize the determinant of the variance-covariance matrix..
//   // This function computes the log determinant of the information matrix using a Choleski decomposition.
//   // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
//   // Input:
//   //     X: design cube of dimensions (q, J, S)
//   //     beta: parameter vector. Must be of length m, with m = (q^3 + 5*q)/6
//   //     verbose: integer that expresses the level of verbosity. Mainly used in other functions and too much useful by itself.
//
//   double eff_crit = 0.0; // We want to minimize this
//
//   int m = beta_mat.n_cols;
//   int n_sims = beta_mat.n_rows;
//
//   arma::mat I(m-1, m-1, fill::zeros);
//
//   // Accumulator
//   double acc = 0.0;
//
//   // Flag in case there's an error in the Cholesky decomposition
//   bool error_flag = false;
//
//   // Iterate over all prior draws
//   for(int i = 0; i < n_sims; i++){
//     I = getInformationMatrixMNL(X, beta); // beta now is a vector. Should put a loop here.
//     if(verbose >= 5) Rcout << "Information matrix. I = \n" << I << std::endl;
//
//
//     arma::mat L;
//     try{
//       L = chol(I);
//       acc = acc - 2*sum(log(L.diag()));
//     }
//     catch(const std::runtime_error& e){
//       // If Cholesky decomposition fails, it is likely because information matrix
//       // was not numerically positive definite.
//       // If this happens, it is probably because a numerical inestability.
//       // The function then returns the efficiency value as a big positive number, this
//       // way the algorithm does nothing in this iteration because the algorithm thinks
//       // there was no improvement when swapping the proportions.
//
//       // If Cholesky decomposition failed, returns a flag to return a final value of 100
//       error_flag = true;
//       Rcpp::warning("Error in Cholesky decomposition with message: ", e.what(), "\nReturning eff_crit = ", eff_crit);
//       break;
//     }
//   }
//
//   if(err_flag){
//     eff_crit = 100;
//   } else{
//     eff_crit = acc/n_sums;
//   }
//
//
//   return eff_crit;
// }
//
//
//
//
//
// // [[Rcpp::export]]
// double getICritValueMNL(arma::cube& X, arma::mat& beta_mat, int verbose, arma::mat& W){
//   int q = X.n_rows;
//
//   int m = beta_mat.n_cols;
//
//   arma::mat I(m-1, m-1, fill::zeros);
//   I = getInformationMatrixMNL(X, beta); // beta now is a vector. Should put a loop here.
//
//   if(verbose >= 5) Rcout << "Information matrix. I = \n" << I << std::endl;
//
//   double eff_crit; // We want to minimize this
//
//
//   // Attempt to do a Cholesky decomposition on the information matrix
//   arma::mat L;
//   try{
//     L = chol(I);
//
//     arma::mat A = solve(trimatl(L.t()), W);
//     arma::mat C = solve(trimatu(L), A);
//     eff_crit = log(trace(C));
//   }
//   catch(const std::runtime_error& e){
//     // If Cholesky decomposition fails, it is likely because information matrix
//     // was not numerically positive definite.
//     // If this happens, it is probably because a numerical inestability.
//     // The function then returns the efficiency value as a big positive number, this
//     // way the algorithm does nothing in this iteration because the algorithm thinks
//     // there was no improvement when swapping the proportions.
//     eff_crit = 100;
//     Rcpp::warning("Error in Cholesky decomposition with message: ", e.what(), "\nReturning eff_crit = ", eff_crit);
//   }
//
//   return eff_crit;
// }





// [[Rcpp::export]]
double getOptCritValueMNL(arma::cube& X, arma::mat& beta_mat, int verbose, int opt_crit, arma::mat& W){
  //  opt_crit: optimality criterion: 0 (D-optimality) or 1 (I-optimality)


  int m = beta_mat.n_cols;
  int n_sims = beta_mat.n_rows;



  // Initialize information matrix with zeros
  arma::mat I(m-1, m-1, fill::zeros);

  // Final efficiency criterion value. We want to minimize this.
  double eff_crit_val = 0.0;

  // Accumulator
  double acc = 0.0;

  // Flag in case there's an error in the Cholesky decomposition
  bool error_flag = false;

  vec beta = zeros(m-1);

  // Iterate over all prior draws
  // This seems to introduce overhead when there is only one row. It takes twice as long as before.
  // I don't really get why it is longer. It's not the for loop or the creation of the vector from the
  // row, since I've tested that.
  // Don't know what's going on.
  for(int i = 0; i < n_sims; i++){
    beta = conv_to<vec>::from(beta_mat.row(i));
    I = getInformationMatrixMNL(X, beta);
    if(verbose >= 5) Rcout << "Information matrix. I = \n" << I << std::endl;


    arma::mat L;
    try{
      L = chol(I);

      if(opt_crit == 0){
        // D-optimality
        acc = acc - 2*sum(log(L.diag()));
      } else{
        // I-optimality
        arma::mat A = solve(trimatl(L.t()), W);
        arma::mat C = solve(trimatu(L), A);
        acc = acc + log(trace(C));
      }
    }

    catch(const std::runtime_error& e){
      // If Cholesky decomposition fails, it is likely because information matrix
      // was not numerically positive definite.
      // If this happens, it is probably because a numerical inestability.
      // The function then returns the efficiency value as a big positive number, this
      // way the algorithm does nothing in this iteration because the algorithm thinks
      // there was no improvement when swapping the proportions.

      // If Cholesky decomposition failed, returns a flag to return a final value of 100
      error_flag = true;
      Rcpp::warning("Error in Cholesky decomposition with message: ", e.what(), "\nReturning eff_crit_val = ", eff_crit_val);
      break;
    }
  }

  if(error_flag){
    eff_crit_val = 100;
  } else{
    eff_crit_val = acc/n_sims;
  }

  return eff_crit_val;
}







// [[Rcpp::export]]
arma::cube findBestCoxDirMNLDiscrete(arma::mat& cox_dir, arma::cube& X_in, arma::mat& beta_mat,
                             int k, int s, double opt_crit_value_best,
                             int verbose, int opt_crit, arma::mat& W) {
  // Function that returns the design that minimizes the optimality criterion value.
  // Returns a cube of dimension (q, J, S) with a design that minimizes the value of the
  //         optimality criterion value.
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal
  //         choice designs for mixtures (2017)
  // Input:
  //     cox_dir: Matrix with Cox direction with q columns. Each row sums up to 1.
  //     X_in: design cube of dimensions (q, J, S).
  //     beta_mat: parameter matrix. If not Bayesian, must have m rows, with m = (q^3 + 5*q)/6.
  //     k: Cox direction index (1 to q).
  //     s: integer s, corresponding to a choice set in 1 to S.
  //     opt_crit_value_best: Efficiency value with which the new efficiencies are compared to.
  //     verbose: integer that expresses the level of verbosity. Mainly used in other functions and
  //           too much useful by itself.
  //     opt_crit: optimality criterion: 0 (D-optimality) or 1 (I-optimality)
  //     W: moment matrix

  // // Create new cube, otherwise it is modified in R too
  arma::cube X = X_in;

  int q = X.n_rows;

  arma::vec x_k(q);
  double opt_crit_value_j;
  int n_cox_points = cox_dir.n_rows;

  for(int j = 0; j < n_cox_points; j++){
    x_k = X(arma::span::all, arma::span(k-1), arma::span(s-1));
    // same as:
    // X(arma::span::all, arma::span(k-1), arma::span(s-1)) = cox_dir.row(j);
    // if the previous worked
    // Basically, here we replace the j-th row of the Cox direction matrix in the corresponding slice and row of X
    for(int l = 0; l < q; l++){
      X(l, k-1, s-1) = cox_dir(j, l);
    }

    if(verbose >= 4){
      Rcout << "\tj = " << j << " (of " << n_cox_points << "), ";
    }

    opt_crit_value_j = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W);

    if(verbose >= 4){
      Rcout << "opt_crit_value_j = "  << opt_crit_value_j << "\n";
    }

    if(verbose >= 6){
      Rcout << "X (with swapped vector) = \n"  << X << "\n";
    }

    //  The D-optimality criterion seeks to maximize the determinant of the information matrix.
    // If new criterion value is better, then keep the new one. If it's not, keep the old one.
    if(opt_crit_value_j < opt_crit_value_best) {
      // This design has a better criterion value, so we keep the design and update the best value
      opt_crit_value_best = opt_crit_value_j;
    } else{
      // This design does not have a better criterion value, so we return to the old design.
      for(int l = 0; l < q; l++){
        X(l, k-1, s-1) = x_k(l);
      }
    }
  }
  return X;
}




// [[Rcpp::export]]
arma::cube changeIngredientDesignMNL(double theta, arma::cube& X, int i, int j, int s){
  // Returns a new design cube Y changing the i-th ingredient in the j-th
  // alternative in s-th choice set is changed to theta in the original design cube X.
  // Since theta is an ingredient proportion, it must be between 0 and 1.
  // Indices i, j and k are 0-indexed.


  // Create new cube Y that is identical to the one pointed by X.
  // Note: This may not be the most computationally effective way, but ut's the easiest at the moment.
  arma::cube Y = X;

  // Number of ingredients
  int q = Y.n_rows;

  // Create a vector with the ingredients of the j-th alternative of the s-th choice set.
  arma::vec x = X(arma::span::all, arma::span(j), arma::span(s));

  double delta = theta - x(i);

  // recompute proportions:
  vec setDiff_aux = linspace<vec>(0, q-1, q);
  vec setDiff = removeElement(setDiff_aux, i);
  int k;
  double result;

  for(int k_aux = 0; k_aux < setDiff.n_elem; k_aux++){
    k = setDiff(k_aux);

    if(abs(1 - x(i)) < 1e-16) {
      // In case x(i) is numerically 1, it will return a numeric zero such that the vector sums up to 1
      // Almost the same as doing result = 0;
      result = (1 - x(i))/(q-1);
    } else{
      // Other case
      result = x(k) - delta*x(k)/(1 - x(i));
    }

    x(k) = result;
  }

  x(i) = theta;

  // Replace the design Y with the recomputed proportions according to Cox direction
  for(int l = 0; l < q; l++){
    Y(l, j, s) = x(l);
  }

  return(Y);

}




// [[Rcpp::export]]
double efficiencyCoxScheffeMNL(double theta, arma::cube& X, arma::mat& beta_mat,
                               int i, int j, int s,
                               int opt_crit, arma::mat& W){
  // Computes efficiency criterion of a design cube X but where the i-th ingredient in the j-th
  // alternative in s-th choice set is changed to theta.
  // Since theta is an ingredient proportion, it must be between 0 and 1.
  // Indices i, j and k are 0-indexed.
  // We want to minimize this.

  arma::cube Y = changeIngredientDesignMNL(theta, X, i, j, s);


  // Return utility function value. We want to minimize this.
  return(getOptCritValueMNL(Y, beta_mat, 0, opt_crit, W));
}



arma::cube findBestCoxDirMNLBrent(
    arma::cube& X, arma::mat& beta_mat, int i, int j, int s, int opt_crit, arma::mat& W,
    double lower = 0, double upper = 1, double tol = 0.0001) {

  auto f = [&X, &beta_mat, i, j, s, opt_crit, &W](double theta){
    return efficiencyCoxScheffeMNL(theta, X, beta_mat, i, j, s, opt_crit, W);
  };

  double theta_brent, theta_star, f_star;
  double f_brent = brent::local_min_mb(lower, upper, tol, f, theta_brent);

  // Check end points:
  double f_lower = efficiencyCoxScheffeMNL(lower, X, beta_mat, i, j, s, opt_crit, W);
  double f_upper = efficiencyCoxScheffeMNL(upper, X, beta_mat, i, j, s, opt_crit, W);

  f_star = std::min({f_lower, f_upper, f_brent});

  if(f_star == f_brent){
    theta_star = theta_brent;
  } else{
    if(f_star == f_lower){
      theta_star = lower;
    } else{
      theta_star = upper;
    }
  }

  return(changeIngredientDesignMNL(theta_star, X, i, j, s));

}











// To implement the Bayesian part, beta will have to be a vector of prior simulations

// [[Rcpp::export]]
Rcpp::List mixtureCoordinateExchangeMNL(
    arma::cube X_orig, arma::mat beta_mat, int max_it, int verbose, int opt_crit, arma::mat W,
    int opt_method, double lower, double upper, double tol, int n_cox_points){
  // Performs the coordinate exchange algorithm for a Multinomial Logit Scheffé model.
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
  // X: 3 dimensional cube with dimensions (q, J, S) where:
  //    q is the number of ingredient proportions
  //    J is the number of alternatives within a choice set
  //    S is the number of choice sets
  // beta_mat: Matrix of parameters. m = (q^3 + 5*q)/6.
  // n_cox_points: Number of points to use in the discretization of Cox direction.
  // max_it: Maximum number of iteration that the coordinate exchange algorithm will do.
  // verbose: level of verbosity. 6 levels, in which level prints the previous plus additional things:
  //    1: Print the optimality criterion value in each iteration and a final summary
  //    2: Print the values of k, s, i, and optimality criterion value in each subiteration
  //    3: Print the resulting X after each iteration, i.e., after each complete pass on the data
  //    4: Print optimality criterion value for each point in the Cox direction discretization
  //    5: Print the resulting X and information matrix after each subiteration
  //    6: Print the resulting X or each point in the Cox direction discretization
  //    opt_crit: optimality criterion: 0 (D-optimality) or 1 (I-optimality)
  //    W: moment matrix
  // Returns an Rcpp::List object with the following objects:
  //    X_orig: The original design. Cube with dimensions (q, J, S).
  //    X: The optimized design. Cube with dimensions (q, J, S).
  //    opt_crit_value_orig: Optimality criterion value of the original design.
  //    opt_crit_value: Optimality criterion value of the optimized design.
  //    n_iter: Number of iterations performed.

  // Does not do input checks because the R wrapper function does them.

  // Create new cube, otherwise it is modified in R too
  arma::cube X = X_orig;

  int J = X.n_cols;
  int q = X.n_rows;
  int S = X.n_elem/(X.n_cols*X.n_rows);

  // Create matrix with appropriate dimensions for Cox direction in each iteration
  arma::mat cox_dir(n_cox_points, q);

  // Vector of ingredient proportions
  arma::vec x(q);

  double opt_crit_value_orig = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W);
  double opt_crit_value_best = opt_crit_value_orig;
  double opt_crit_value_aux = 1e308; // +Inf

  // Coordinate exchanges
  int it = 0;
  while(it < max_it){

    // Checking interruption every 2 iterations
    // As in https://teuder.github.io/rcpp4everyone_en/270_miscellaneous.html
    if (it % 2 == 0){
      Rcpp::checkUserInterrupt();
    }

    it = it + 1;
    if(verbose >= 1) Rcout << "Iter: " << it << ", Optimality criterion value: " << opt_crit_value_best << std::endl;

    // If there was no improvement in this iteration
    if(abs(opt_crit_value_aux - opt_crit_value_best) < 1e-16) break;

    opt_crit_value_aux = opt_crit_value_best;

    for(int k = 1; k <= J; k++){
      // if(verbose >= 2) Rcout << "k = " << k << std::endl;

      for(int s = 1; s <= S; s++){
        // if(verbose >= 2) Rcout << "\ts = " << s << std::endl;

        for(int i = 0; i < q; i++){
          // if(verbose >= 2) Rcout << "\t\ti = " << i << std::endl;
          if(verbose >= 2) Rcout << "\nIter: " << it <<  ", k = " << k << ", s = " << s << ", i = " << i << std::endl;

          // populate x vector with corresponding ingredient proportions
          for(int l = 0; l < q; l++){
            x(l) = X(l, k-1, s-1);
          }

          if(opt_method == 0){
            X = findBestCoxDirMNLBrent(X, beta_mat, i, k-1, s-1, opt_crit, W, lower, upper, tol);
          } else{
            cox_dir = computeCoxDirection(x, i+1, n_cox_points, verbose);
            X = findBestCoxDirMNLDiscrete(cox_dir, X, beta_mat, k, s, opt_crit_value_best, verbose, opt_crit, W);
          }

          opt_crit_value_best = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W);

          if(verbose >= 2) Rcout << "Opt-crit-value: " << opt_crit_value_best << std::endl;

          if(verbose >= 5){
            Rcout << "X =\n" << X << std::endl;
          }

        } // end for i

      } // end for s

    } // end for k

    if(verbose >= 3) Rcout << "X =\n" << X << std::endl;

    if(verbose >= 2) Rcout << std::endl << std::endl;

  } // end while

  if(verbose >= 1){
    Rcout << std::endl;
    Rcout << "Original Optimality criterion value: " << opt_crit_value_orig;
    Rcout << std::endl;
    Rcout << "Final Optimality criterion value: " << opt_crit_value_best;
    Rcout << std::endl;
    Rcout << "Number of iterations: " << it;
    Rcout << std::endl;
  }

  // return object
  return Rcpp::List::create(
    _["X_orig"] = X_orig,
    _["X"] = X,
    _["opt_crit_value_orig"] = opt_crit_value_orig,
    _["opt_crit_value"] = opt_crit_value_best,
    _["n_iter"] = it
  );

} // end function
















// // [[Rcpp::export]]
// List BrentCoxScheffeGaussianNoEdge(arma::mat& X, int j, int i, int order, int opt_crit, arma::mat& W,
//                              double lower = 0, double upper = 1, double tol = 0.0001){
//   auto f = [&X, j, i, order, opt_crit, &W](double theta){
//     return efficiencyCoxScheffeGaussian(theta, X, j, i, order, opt_crit, W);
//     // double efficiencyCoxScheffeGaussian(double theta, arma::mat& X, int j, int i, int order, int opt_crit, arma::mat& W)
//   };
//
//   double theta_star;
//   double f_star = brent::local_min_mb(lower, upper, tol, f, theta_star);
//   return List::create(Named("minimizer") = theta_star, Named("objective_func") = f_star);
// }
//
//
//
//
//
// // [[Rcpp::export]]
// List BrentCoxScheffeGaussian(arma::mat& X, int j, int i, int order, int opt_crit, arma::mat& W,
//                              double lower = 0.0, double upper = 1.0, double tol = 0.0001){
//   auto f = [&X, j, i, order, opt_crit, &W](double theta){
//     return efficiencyCoxScheffeGaussian(theta, X, j, i, order, opt_crit, W);
//   };
//
//   double theta_brent, theta_star, f_star;
//   double f_brent = brent::local_min_mb(lower, upper, tol, f, theta_brent);
//   // Check edge cases
//   double f_lower = efficiencyCoxScheffeGaussian(lower, X, j, i, order, opt_crit, W);
//   double f_upper = efficiencyCoxScheffeGaussian(upper, X, j, i, order, opt_crit, W);
//
//   f_star = std::min({f_lower, f_upper, f_brent});
//
//   if(f_star == f_brent){
//     theta_star = theta_brent;
//   } else{
//     if(f_star == f_lower){
//       theta_star = lower;
//     } else{
//       theta_star = upper;
//     }
//   }
//   return List::create(Named("minimizer") = theta_star, Named("objective_func") = f_star);
// }
//
//
//
//
//
// // [[Rcpp::export]]
// List BrentGloCoxScheffeGaussian(
//     arma::mat& X, int j, int i, int order, int opt_crit, arma::mat& W,
//     double lower = 0, double upper = 1,
//     double initial_guess = 0.5,
//     double hessian_bound = 1e5,
//     double abs_err_tol = 0.0001,
//     double tol = 0.0001){
//   auto f = [&X, j, i, order, opt_crit, &W](double theta){
//     return efficiencyCoxScheffeGaussian(theta, X, j, i, order, opt_crit, W);
//     // double efficiencyCoxScheffeGaussian(double theta, arma::mat& X, int j, int i, int order, int opt_crit, arma::mat& W)
//   };
//
//   double theta_star;
//   double f_star = brent::glomin_mb(
//     lower, upper, initial_guess, hessian_bound, abs_err_tol, tol, f, theta_star);
//   return List::create(Named("minimizer") = theta_star, Named("objective_func") = f_star);
// }










// // Banana function
//
// // [[Rcpp::export]]
// double banana_xy(double x, double y){
//   return (1 - x)*(1-x) + 100*(y - x*x)*(y - x*x);
// }
//
//
// std::function<double(double)> create_banana(double y) {
//   auto out = [y](double x) {
//     return(banana_xy(x, y));
//   };
//   return out;
// }
//
// // [[Rcpp::export]]
// double banana_x_y1(double x){
//   std::function<double(double)> foo = create_banana(1.0);
//   return(foo(x));
// }
//
// // [[Rcpp::export]]
// double banana_x_y2(double x){
//   std::function<double(double)> foo = create_banana(2.0);
//   return(foo(x));
// }
//
//
//
// // [[Rcpp::export]]
// List min_banana_x_y1(double lower = -10, double upper = 10, double tol = 0.0001){
//   double x_star;
//   double f_star = brent::local_min_mb(lower, upper, tol, banana_x_y1, x_star);
//   return List::create(Named("minimizer") = x_star, Named("objective_func") = f_star);
// }
//
//
//
//
// // [[Rcpp::export]]
// List minimize_banana_fixed_y(double y = 1.0, double lower = -10, double upper = 10, double tol = 0.0001){
//   auto f = [y](double x){ return banana_xy(x, y); };
//
//   double x_star;
//   double f_star = brent::local_min_mb(lower, upper, tol, f, x_star);
//   return List::create(Named("minimizer") = x_star, Named("objective_func") = f_star);
// }
//
//
//
//
//
// // [[Rcpp::export]]
// double banana_xy2(double x, double &y){
//   return (1 - x)*(1-x) + 100*(y - x*x)*(y - x*x);
// }
//
//
// // [[Rcpp::export]]
// List minimize_banana_fixed_y2(double &y,
//                               double lower = -10, double upper = 10, double tol = 0.0001){
//   auto f = [&y](double x){ return banana_xy2(x, y); };
//
//   double x_star;
//   double f_star = brent::local_min_mb(lower, upper, tol, f, x_star);
//   return List::create(Named("minimizer") = x_star, Named("objective_func") = f_star);
// }








