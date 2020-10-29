#include <RcppArmadillo.h>
#include <functional>
#include "utils.hpp"
#include "brent.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;


////////////////////////////////////////
// Multinomial logit (MNL) model
////////////////////////////////////////



// [[Rcpp::export]]
arma::mat getXsMNL(arma::cube& X, int s, int order){
  // Function that returns the design matrix of choice set s.
  // Final matrix is of dimension (J, m-1) with m = (q^3 + 5*q)/6 (if 3rd ordder Scheffe)
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
  // Input should be
  //     X: design cube of dimensions (q, J, S) integer s, corresponding to a choice set in 1 to S.
  //     order: order of Scheffe model
  int q = X.n_rows;
  int J = X.n_cols;
  int m;


  if(order == 1){
    m = q;
  } else{
    if(order == 2){
      m = q*(q-1)/2 + q;
    } else{
      m = (q*q*q + 5*q)/6; // = q + q*(q-1)/2 + q*(q-1)*(q-2)/6 = (q^3+ 5*q)/6
    }
  }

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

  if(order >= 2){
    // second part
    for(int i = 0; i < q-1; i++){
      for(int k = i+1; k < q; k++){
        for(int j = 0; j < J; j++){
          Xs(j, col_counter) = X(i, j, s-1)*X(k, j, s-1);
        }
        col_counter++;
      }
    }
  }


  if(order >= 3){
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
  // if(m != beta.n_elem) stop("Incompatible q in beta and X");

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
arma::mat getInformationMatrixMNL(arma::cube& X, arma::vec& beta, int order){
  // Function that returns the information matrix for design cube X and parameter vector beta.
  // It is the sum of the information matrices of the S choice sets.
  // Final matrix is of dimension (m-1, m-1), with m = (q^3 + 5*q)/6
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
  // Input:
  //     X: design cube of dimensions (q, J, S)
  //     beta: parameter vector. Must be of length m, with m = (q^3 + 5*q)/6 (if 3rd order Scheffe)
  //     order: order of Scheffe model

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
    Xs = getXsMNL(X, s, order);
    ps = getPsMNL(X, beta, s, Xs);;
    ps_ps_t = ps*ps.t();
    middle = ps_ps_t;
    middle.diag() = ps_ps_t.diag() - ps;
    I = I - (Xs.t())*middle*Xs;
  }
  return I;
}





// [[Rcpp::export]]
double getOptCritValueMNL(arma::cube& X, arma::mat& beta_mat, int verbose, int opt_crit, arma::mat& W, int order){
  //  opt_crit: optimality criterion: 0 (D-optimality) or 1 (I-optimality)


  int m = beta_mat.n_cols;
  int n_sims = beta_mat.n_rows;



  // Initialize information matrix with zeros
  arma::mat I(m-1, m-1, fill::zeros);

  // Final efficiency criterion value. We want to minimize this.
  double eff_crit_val = 0.0;

  // Accumulator
  double acc = 0.0;

  vec beta = zeros(m-1);

  arma::mat L; // matrix for Cholesky decomposition
  bool chol_success; // flag for success of the decomposition

  // Iterate over all prior draws
  for(int i = 0; i < n_sims; i++){
    beta = conv_to<vec>::from(beta_mat.row(i));
    I = getInformationMatrixMNL(X, beta, order);
    if(verbose >= 5) Rcout << "Information matrix. I = \n" << I << std::endl;

    // If the decomposition fails chol(R,X) resets R and returns a bool set to false (exception is not thrown) (http://arma.sourceforge.net/docs.html#chol)
    chol_success = chol(L, I);

    if(chol_success){
      if(opt_crit == 0){
        // D-optimality
        acc = acc - 2*sum(log(L.diag()));
      } else{
        // I-optimality
        arma::mat A = solve(trimatl(L.t()), W);
        arma::mat C = solve(trimatu(L), A);
        acc = acc + log(trace(C));
      }
    } else {
      // If Cholesky decomposition fails, it is likely because information matrix
      // was not numerically positive definite.
      // If this happens, it is probably because a numerical inestability.
      // The function then returns the efficiency value as a big positive number, this
      // way the algorithm does nothing in this iteration because the algorithm thinks
      // there was no improvement when swapping the proportions.

      // If Cholesky decomposition failed, exits and returns a final value of 1000
      // Rcpp::warning("Error in Cholesky decomposition with message: '%c'.\n\tReturning eff_crit = 1000.");
      break;
    }
  } // end for

  if(!chol_success){
    eff_crit_val = 1000;
  } else{
    eff_crit_val = acc/n_sims;
  }

  return eff_crit_val;
}







// [[Rcpp::export]]
void findBestCoxDirMNLDiscrete( // void but modifies X
    arma::mat& cox_dir, arma::cube& X, arma::mat& beta_mat,
    int k, int s, double opt_crit_value_best,
    int verbose, int opt_crit, arma::mat& W, int order) {
  // Function that returns the design that minimizes the optimality criterion value.
  // Returns a cube of dimension (q, J, S) with a design that minimizes the value of the
  //         optimality criterion value.
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal
  //         choice designs for mixtures (2017)
  // Input:
  //     cox_dir: Matrix with Cox direction with q columns. Each row sums up to 1.
  //     X_in: design cube of dimensions (q, J, S).
  //     beta_mat: parameter matrix. If not Bayesian, must have m rows, with m = (q^3 + 5*q)/6. (if order 3 Scheffe)
  //     k: Cox direction index (1 to q).
  //     s: integer s, corresponding to a choice set in 1 to S.
  //     opt_crit_value_best: Efficiency value with which the new efficiencies are compared to.
  //     verbose: integer that expresses the level of verbosity. Mainly used in other functions and
  //           too much useful by itself.
  //     opt_crit: optimality criterion: 0 (D-optimality) or 1 (I-optimality)
  //     W: moment matrix
  //     order: order of Scheffe model

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

    opt_crit_value_j = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W, order);

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
}





// [[Rcpp::export]]
void changeIngredientDesignMNL(double theta, arma::cube& X, int i, int j, int s){ // Void but modifies X
  // Modifies cube X changing the i-th ingredient in the j-th alternative in s-th choice set to theta.
  // Since theta is an ingredient proportion, it must be between 0 and 1.
  // Indices i, j and k are 0-indexed.


  // Number of ingredients
  double q = X.n_rows;

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

    if(abs(1 - x(i)) < 1e-16) { // In case x(i) is numerically 1
      if(abs(1 - x(i)) < 1e-16) { // In case x(i) is numerically 1
        if(abs(delta + 1) < 1e-16){
          // If delta is numerically -1, it means that the change is from 1 to 0.
          // Then, the rest of the ingredients must be 1/(q-1)
          result = 1.0/(q - 1.0);
        } else{
          // If delta is not -1 and x(i) is numerically 1, then
          // it means that the rest of the ingredients were 0 and it should remain that way.
          // Almost the same as doing result = 0, but this is safer numerically.
          result = (1.0 - x(i))/(q - 1.0);
        }
      } else{ // In case x(i) is not numerically 1
        // Other case
        result = x(k) - delta*x(k)/(1 - x(i));
      }
    } else{ // In case x(i) is not numerically 1
      result = x(k) - delta*x(k)/(1 - x(i));
    }
    x(k) = result;
  }

  x(i) = theta;

  // Replace the design X with the recomputed proportions according to Cox direction
  for(int l = 0; l < q; l++){
    X(l, j, s) = x(l);
  }

}




// [[Rcpp::export]]
double efficiencyCoxScheffeMNL(double theta, arma::cube& X, arma::mat& beta_mat,
                               int i, int j, int s,
                               int opt_crit, arma::mat& W, int order){
  // Computes efficiency criterion of a design cube X but where the i-th ingredient in the j-th
  // alternative in s-th choice set is changed to theta.
  // Since theta is an ingredient proportion, it must be between 0 and 1.
  // Indices i, j and k are 0-indexed.
  // We want to minimize this.

  double q = X.n_rows;

  // Temporarily store the corresponding row in a vector
  arma::vec x_row(q);
  for(int l = 0; l < q; l++){
    x_row(l) = X(l, j, s);
  }

  changeIngredientDesignMNL(theta, X, i, j, s);

  // Utility function value. We want to minimize this.
  double utility_funct_value = getOptCritValueMNL(X, beta_mat, 0, opt_crit, W, order);

  // return X to its original value
  for(int l = 0; l < q; l++){
    X(l, j, s) = x_row(l);
  }

  return(utility_funct_value);

}



void findBestCoxDirMNLBrent( // void but modifies X
    arma::cube& X, arma::mat& beta_mat, int i, int j, int s, int opt_crit, int order, arma::mat& W,
    double lower = 0, double upper = 1, double tol = 0.0001) {
  // Finds the best point in the Cox direction using Brent's optimization method


  // The optimality criterion value in the design that is being used as input
  double f_original = getOptCritValueMNL(X, beta_mat, 0, opt_crit, W, order);


  // Helper function that depends only on theta (the ingredient proportion that is being changed now)
  auto f = [&X, &beta_mat, i, j, s, opt_crit, &W, order](double theta){
    return efficiencyCoxScheffeMNL(theta, X, beta_mat, i, j, s, opt_crit, W, order);
  };

  // theta_brent: the "optimal" theta that is going to be returned by Brent's method
  // theta_star: the actual optimal theta
  // f_star: the value of the objective function evaluated in the optimal theta
  double theta_brent, theta_star, f_star;
  double f_brent = brent::local_min_mb(lower, upper, tol, f, theta_brent);

  // Check end points:
  double f_lower = efficiencyCoxScheffeMNL(lower, X, beta_mat, i, j, s, opt_crit, W, order);
  double f_upper = efficiencyCoxScheffeMNL(upper, X, beta_mat, i, j, s, opt_crit, W, order);

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

  // If Brent's method didn't do any improvement, leave the original design as it is.
  // If it did, replace the design with theta star
  if(f_original > f_star){
    changeIngredientDesignMNL(theta_star, X, i, j, s);
  }
}












// [[Rcpp::export]]
Rcpp::List mixtureCoordinateExchangeMNL(
    arma::cube X_orig, arma::mat beta_mat, int order, int max_it, int verbose, int opt_crit, arma::mat W,
    int opt_method, double lower, double upper, double tol, int n_cox_points){
  // See mnl_mixture_coord_exch() in R for details.
  // X_orig: If an initial design is to be supplied, thenit must be a 3 dimensional array with dimensions (q, J, S), with q, J, and S are defined above.
  // beta_mat: Prior parameters. For a locally optimal design, it should be a numeric vector of length m = (q^3 + 5*q)/6. For a pseudo-Bayesian design,
  //           it must be a matrix with prior simulations of size (nxm) where m is previously defined and m is the number of prior draws,
  //           i.e., there is a prior draw per row.
  // order: order of Scheffe's model. Must be 1, 2, or 3.
  // max_it: integer for maximum number of iterations that the coordinate exchange algorithm will do
  // verbose: level of verbosity.
  // opt_crit: optimality criterion: D-optimality (0) or I-optimality (1).
  // W: moments matrix
  // opt_method: Optimization method in each step of the coordinate exchange algorithm.
  //             0 for Brent, 1 for Cox's direction discretization.
  // lower: lower bound in Brent's optimization method
  // upper: upper bound in Brent's optimization method
  // tol: A positive error tolerance in Brent's method.
  // n_cox_points: number of points to use in the discretization of Cox direction. Ignored if opt_method is Brent.


  // Does not do input checks because the R wrapper function does them.

  // Create a vector to store the values of the efficiency metric in each iteration.
  arma::vec efficiency_value_per_iteration(max_it + 1, fill::zeros);

  // Create new cube, otherwise it is modified in R too
  arma::cube X = X_orig;

  int J = X.n_cols;
  int q = X.n_rows;
  int S = X.n_elem/(X.n_cols*X.n_rows);

  // Create matrix with appropriate dimensions for Cox direction in each iteration
  arma::mat cox_dir(n_cox_points, q);

  // Vector of ingredient proportions
  arma::vec x(q);

  double opt_crit_value_orig = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W, order);
  double opt_crit_value_best = opt_crit_value_orig;
  double opt_crit_value_aux = 1e308; // +Inf

  // Original efficiency criterion value
  efficiency_value_per_iteration(0) = opt_crit_value_best;

  // Coordinate exchanges
  int it = 0;
  while(it < max_it){

    // Checking interruption every 2 iterations
    // As in https://teuder.github.io/rcpp4everyone_en/270_miscellaneous.html
    if (it % 2 == 0){
      Rcpp::checkUserInterrupt();
    }

    if(verbose >= 1) Rcout << "Iter " << it << ". Optimality criterion value: " << opt_crit_value_best << std::endl;

    // If there was no improvement in this iteration
    if(abs(opt_crit_value_aux - opt_crit_value_best) < 1e-16) break;

    it = it + 1;

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
            findBestCoxDirMNLBrent(X, beta_mat, i, k-1, s-1, opt_crit, order, W, lower, upper, tol);
          } else{
            cox_dir = computeCoxDirection(x, i+1, n_cox_points, verbose);
            findBestCoxDirMNLDiscrete(cox_dir, X, beta_mat, k, s, opt_crit_value_best, verbose, opt_crit, W, order);
          }

          opt_crit_value_best = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W, order);

          if(verbose >= 2) Rcout << "Opt-crit-value: " << opt_crit_value_best << std::endl;

          if(verbose >= 5){
            Rcout << "X =\n" << X << std::endl;
          }

        } // end for i

      } // end for s

    } // end for k

    if(verbose >= 3) Rcout << "X =\n" << X << std::endl;

    if(verbose >= 2) Rcout << std::endl << std::endl;

    efficiency_value_per_iteration(it) = opt_crit_value_best;

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















