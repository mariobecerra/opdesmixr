#include <RcppArmadillo.h>
#include <functional>
#include "utils.hpp"
#include "brent.hpp"


using namespace Rcpp;
using namespace arma;
using namespace std;

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
  // arma::mat X_ord_2 = getScheffeGaussianOrder2(X);

  int n_col_X_ord_2 = q*(q-1)/2 + q;
  int n_col_X_m = (q*q*q + 5*q)/6;

  arma::mat X_m(n, n_col_X_m);

  // Copy X matrix into first q columns of X_m (first order part)
  for(int i = 0; i < n; i++){
    for(int j = 0; j < q; j++){
      X_m(i,j) = X(i,j);
    }
  }

  // Fill the part of the second order
  int k = q-1;
  for(int i = 0; i < (q-1); i++){
    for(int j = i+1; j < q; j++){
      k = k+1;
      for(int row = 0; row < n; row++){
        X_m(row,k) = X(row,i)*X(row,j);
      }
    }
  }

  // Fill the part of the third order
  int l = n_col_X_ord_2 - 1;
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




arma::mat getScheffeGaussianPV22(arma::mat& X, int m = 0){
  // m = number of process variables
  // X = design matrix

  // Equation 6.12 in Peter Goos's and Peter Jone's book "Optimal Design of Experiments, a Case Study Approach"
  // Second order Scheffé model and quadratic effect of process variables (hence the name PV22)


  int q = X.n_cols - m; // q = Number of ingredient proportions
  int n = X.n_rows;
  int n_col_X_m = q + (q-1)*q/2 + q*m + m*(m-1)/2 + m;
  arma::mat X_m(n, n_col_X_m);

  // Copy X matrix into first q columns of X_m
  for(int i = 0; i < n; i++){
    for(int j = 0; j < q; j++){
      X_m(i, j) = X(i, j);
    }
  }

  int count = q-1;
  // Second part of the sum
  for(int i = 0; i < (q-1); i++){
    for(int j = i+1; j < q; j++){
      count++;
      for(int row = 0; row < n; row++){
        X_m(row, count) = X(row, i)*X(row, j);
      }
    }
  }

  // Third part of the sum
  for(int i = 0; i < m; i++){
    for(int k = 0; k < q; k++){
      count++;
      for(int row = 0; row < n; row++){
        X_m(row, count) = X(row, k)*X(row, i+q); //x_k*z_i
      }
    }
  }

  // Fourth part of the sum
  for(int i = 0; i < m-1; i++){
    for(int j = i+1; j < m; j++){
      count++;
      for(int row = 0; row < n; row++){
        X_m(row, count) = X(row, i+q)*X(row, j+q); //z_i*z_j
      }
    }
  }

  // Fifth part of the sum
  for(int i = 0; i < m; i++){
    count++;
    for(int row = 0; row < n; row++){
      X_m(row, count) = X(row, i+q)*X(row, i+q); //z_i^2
    }
  }

  return X_m;
}




// [[Rcpp::export]]
arma::mat getScheffeGaussian(arma::mat& X, int order, int m = 0){

  if(order == 1) return(X);
  if(order == 2) return(getScheffeGaussianOrder2(X));
  if(order == 3) return(getScheffeGaussianOrder3(X));
  if(order == 4) return(getScheffeGaussianPV22(X, m));
  else stop("Error.");

}



// [[Rcpp::export]]
double getDCritValueGaussian(arma::mat& X, int order, int m = 0){
  arma::mat X_m = getScheffeGaussian(X, order, m);
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
double getICritValueGaussian(arma::mat& X, int order, int q, arma::mat& W, int m = 0){

  // if(order != 3) stop("Only special cubic models are allowed (i.e., order = 3).");

  arma::mat X_m = getScheffeGaussian(X, order, m);
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
double getOptCritValueGaussian(arma::mat& X, int order, int q, int opt_crit, arma::mat& W, int m = 0){
  //  opt_crit: optimality criterion: 0 (D-optimality) or 1 (I-optimality)
  if(opt_crit == 0){
    return(getDCritValueGaussian(X, order, m));
  } else{
    return(getICritValueGaussian(X, order, q, W, m));
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







arma::mat changeIngredientDesignCoxGaussian(double theta, arma::mat& X, int i, int j){
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
  double q = Y.n_cols;

  // Create a vector with the i-th row of the design matrix X
  arma::vec x_row(q);
  for(int col = 0; col < q; col++){
    x_row(col) = X(i, col);
  }

  double delta = theta - x_row(j);

  // recompute proportions:
  vec setDiff_aux = linspace<vec>(0, q-1, q);
  vec setDiff = removeElement(setDiff_aux, j);
  int k;

  for(int k_aux = 0; k_aux < setDiff.n_elem; k_aux++){
    k = setDiff(k_aux);

    if(abs(1 - x_row(j)) < 1e-16) { // In case x_row(j) is numerically 1
      if(abs(delta + 1) < 1e-16){
        // If delta is numerically -1, it means that the change is from 1 to 0.
        // Then, the rest of the ingredients must be 1/(q-1)
        x_row(k) = 1/(q - 1.0);
      } else{
        // If delta is not -1 and x_row(j) is numerically 1, then
        // it means that the rest of the ingredients were 0 and it should remain that way.
        // Almost the same as doing x_row(k) = 0, but this is safer numerically.
        x_row(k) = (1 - x_row(j))/(q-1.0);
      }
    } else{ // In case x_row(j) is NOT numerically 1
      x_row(k) = x_row(k) - delta*x_row(k)/(1 - x_row(j));
    }
  }

  x_row(j) = theta;

  // replace x_row row with the recomputed proportions according to Cox direction
  for(int col = 0; col < q; col++){
    Y(i, col) = x_row(col);
  }

  return(Y);

}




// [[Rcpp::export]]
double efficiencyCoxScheffeGaussian(double theta, arma::mat& X, int i, int j, int order,
                                    int opt_crit, arma::mat& W){
  // Computes efficiency criterion of a design matrix X but where the j-th ingredient in the
  // i-th observation is changed to theta.
  // Since theta is an ingredient proportion, it must be between 0 and 1.
  // INdices j and i are 0-indexed.
  // We want to minimize this.


  arma::mat Y = changeIngredientDesignCoxGaussian(theta, X, i, j);
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

  // The optimality criterion value in the design that is being used as input
  double f_original = getOptCritValueGaussian(X, order, X.n_cols, opt_crit, W);

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


  // If Brent's method didn't do any improvement, then return the original design
  if(f_original <= f_star){
    return(X);
  } else{
    return(changeIngredientDesignCoxGaussian(theta_star, X, i, j));
  }

}











double efficiencyPVScheffeGaussian(double theta, arma::mat& X, int i, int j, int order,
                                   int opt_crit, arma::mat& W, arma::mat& Y){
  // Computes efficiency criterion of a design matrix X but where the j-th PV in the
  // i-th observation is changed to theta.
  // Indices j and i are 0-indexed.
  // We want to minimize this.

  int q = X.n_cols;
  Y(i, j) = theta;

  // Return utility function value. We want to minimize this.
  return(getOptCritValueGaussian(Y, order, q, opt_crit, W));
}




void findBestPVGaussianBrent( // void but modifies Y
    arma::mat& X, int i, int j, int order, int opt_crit, arma::mat& W,
    double lower, double upper, double tol, arma::mat& Y) {

  // The optimality criterion value in the design that is being used as input
  double f_original = getOptCritValueGaussian(X, order, X.n_cols, opt_crit, W);

  auto f = [&X, i, j, order, opt_crit, &W, &Y](double theta){
    return efficiencyPVScheffeGaussian(theta, X, i, j, order, opt_crit, W, Y);
  };

  double theta_brent, theta_star, f_star;
  double f_brent = brent::local_min_mb(lower, upper, tol, f, theta_brent);
  // Check end points
  double f_lower = efficiencyPVScheffeGaussian(lower, X, i, j, order, opt_crit, W, Y);
  double f_upper = efficiencyPVScheffeGaussian(upper, X, i, j, order, opt_crit, W, Y);

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


  // If Brent's method didn't do any improvement, then return to the original design
  if(f_original <= f_star){
    Y(i, j) = X(i, j);
  } else{
    Y(i, j) = theta_star;
  }

}






// [[Rcpp::export]]
Rcpp::List mixtureCoordinateExchangeGaussian(
    arma::mat X_orig, int order, int max_it, int verbose, int opt_crit,
    arma::mat W, int opt_method, double lower, double upper, double tol, int n_cox_points){
  // Performs the coordinate exchange algorithm for a Scheffé model with Gaussian errors.
  int m = 0;

  // n_runs: number of runs
  // q: number of ingredient proportions
  // n_random_starts: number or random starts. Defaults to 100.
  // X: User supplied design matrix.
  // order: Order of the Scheffé model (1, 2, or 3).
  // opt_method: Optimization method in each step of the coordinate exchange algorithm.
  //    0 for Brent, 1 for Cox's direction discretization.
  // max_it: integer for maximum number of iterations that the coordinate exchange algorithm will do
  // tol: A positive error tolerance in Brent's method.
  // n_cox_points: number of points to use in the discretization of Cox direction
  // verbose level of verbosity.
  // opt_crit optimality criterion: D-optimality (0) or I-optimality (1)
  // W: moment matrix


  // Does not do input checks because the R wrapper function does them.

  // Create a vector to store the values of the efficiency metric in each iteration.
  arma::vec efficiency_value_per_iteration(max_it + 1, fill::zeros);

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
  double opt_crit_value_aux = 1e308; // +Inf

  // Original efficiency criterion value
  efficiency_value_per_iteration(0) = opt_crit_value_best;

  // Coordinate exchanges
  int it = 0;
  while(it < max_it){

    if(verbose >= 1) Rcout << "Iter " << it << ". Optimality criterion value: " << opt_crit_value_best << std::endl;

    // If there was no improvement in this iteration
    if(abs(opt_crit_value_aux - opt_crit_value_best) < 1e-16) break;

    it = it + 1;

    opt_crit_value_aux = opt_crit_value_best;

    for(int run = 1; run <= n_runs; run++){

      // Checking interruption every 2 iterations
      // As in https://teuder.github.io/rcpp4everyone_en/270_miscellaneous.html
      if (run % 2 == 0){
        Rcpp::checkUserInterrupt();
      }


      // Mixture variables
      for(int ing = 0; ing < q; ing++){

        if(verbose >= 2) Rcout << "\nIter: " << it <<  ", run = " << run << ", ing = " << ing << std::endl;

        // populate x vector with corresponding ingredient proportions
        for(int l = 0; l < q; l++){
          x(l) = X(run-1, l);
        }

        if(opt_method == 0){
          X = findBestCoxDirGaussianBrent(X, run-1, ing, order, opt_crit, W, lower, upper, tol);

        } else{
          cox_dir = computeCoxDirection(x, ing+1, n_cox_points, verbose);
          X = findBestCoxDirGaussianDiscrete(cox_dir, X, run, order, opt_crit_value_best, opt_crit, W);
        }


        opt_crit_value_best = getOptCritValueGaussian(X, order, q, opt_crit, W);

        if(verbose >= 2) Rcout << "Opt-crit-value: " << opt_crit_value_best << std::endl;

        if(verbose >= 5){
          Rcout << "X =\n" << X << std::endl;
        }

      } // end for ing

      // Process variables. To do.
      if(m > 0){
        arma::mat Y = X;
        for(int pv = q; pv < m + q; pv++){
          findBestPVGaussianBrent(X, run-1, pv, order, opt_crit, W, lower, upper, tol, Y);
          X(run, pv) = Y(run, pv);
        }
      }



    } // end for run

    efficiency_value_per_iteration(it) = opt_crit_value_best;

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
    _["efficiency_value_per_iteration"] = efficiency_value_per_iteration.head(it + 1)
  );

} // end function


