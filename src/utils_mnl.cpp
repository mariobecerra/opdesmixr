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
arma::mat getXsMNL(arma::cube& X, int s, int order, int n_pv = 0, bool no_choice = false){
  // Function that returns the design matrix of choice set s.
  // For more information about this function, see R's wrapper function mnl_get_Xs()
  int q = X.n_rows - n_pv;
  int J = X.n_cols;
  int m = 0; // initialize to silence warning


  if(order == 1){
    m = q;
  } else{
    if(order == 2){
      m = q*(q-1)/2 + q;
    } else{
      if(order == 3){
        m = (q*q*q + 5*q)/6; // = q + q*(q-1)/2 + q*(q-1)*(q-2)/6 = (q^3+ 5*q)/6
      } else{
        if(order == 4){
          m = q + (q-1)*q/2 + q*n_pv + n_pv*(n_pv-1)/2 + n_pv;
        } else{
          // Error. Not necessary to check here because there are input checks in R.
          stop("Wrong order.");
        }
      }
    }
  }

  int nrow_Xs, ncol_Xs;
  if(!no_choice){
    nrow_Xs = J;
    ncol_Xs = m-1;
  } else{
    nrow_Xs = J+1;
    ncol_Xs = m;
  }

  // Initialize array with zeros
  arma::mat Xs(nrow_Xs, ncol_Xs, fill::zeros);

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


  if(order == 3){
    // third part (no process variables)
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



  if(order == 4){
    // Third part of the sum (with process variables)
    for(int i = 0; i < n_pv; i++){
      for(int k = 0; k < q; k++){

        for(int j = 0; j < J; j++){
          Xs(j, col_counter) = X(k, j, s-1)*X(i+q, j, s-1); //x_k*z_i
        }
        col_counter++;
      }
    }

    // Fourth part of the sum (with process variables)
    for(int i = 0; i < n_pv-1; i++){
      for(int k = i+1; k < n_pv; k++){

        for(int j = 0; j < J; j++){
          Xs(j, col_counter) = X(i+q, j, s-1)*X(k+q, j, s-1); //z_i*z_k
        }
        col_counter++;
      }
    }

    // Fifth part of the sum (with process variables)
    for(int i = 0; i < n_pv; i++){

      for(int j = 0; j < J; j++){
        Xs(j, col_counter) = X(i+q, j, s-1)*X(i+q, j, s-1); //z_i^2
      }
      col_counter++;
    }
  }


  // Fill the 1 for the no choice option
  if(no_choice) Xs(J, m-1) = 1;



  return Xs;
}




// // [[Rcpp::export]]
// arma::mat getXsMNL2(arma::cube& X, int s, int order, int n_pv = 0){
//
//   // TEst to see if including the q-th component changed the design values compared with JMP. It didn't.
//
//   // Function that returns the design matrix of choice set s.
//   // Final matrix is of dimension (J, m-1) with m = (q^3 + 5*q)/6 (if 3rd ordder Scheffe)
//   // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
//   // Input should be
//   //     X: design cube of dimensions (q, J, S) integer s, corresponding to a choice set in 1 to S.
//   //     order: order of Scheffe model
//   int q = X.n_rows - n_pv;
//   int J = X.n_cols;
//   int m = 0; // initialize to silence warning
//
//
//   if(order == 1){
//     m = q;
//   } else{
//     if(order == 2){
//       m = q*(q-1)/2 + q;
//     } else{
//       if(order == 3){
//         m = (q*q*q + 5*q)/6; // = q + q*(q-1)/2 + q*(q-1)*(q-2)/6 = (q^3+ 5*q)/6
//       } else{
//         if(order == 4){
//           m = q + (q-1)*q/2 + q*n_pv + n_pv*(n_pv-1)/2 + n_pv;
//         } else{
//           // Error. Not necessary to check here because there are input checks in R.
//           stop("Wrong order.");
//         }
//       }
//     }
//   }
//
//   // Initialize array with zeros
//   arma::mat Xs(J, m, fill::zeros);
//
//   // Column counter. Equation has three terms, so the final matrix is populated in three parts.
//   int col_counter = 0;
//
//   // First part
//   for(int i = 0; i < q; i++){
//     for(int j = 0; j < J; j++){
//       // subtract 1 from s because it is 1-indexed
//       Xs(j, col_counter) = X(i, j, s-1);
//     }
//     col_counter++;
//   }
//
//   if(order >= 2){
//     // second part
//     for(int i = 0; i < q; i++){
//       for(int k = i+1; k < q; k++){
//         for(int j = 0; j < J; j++){
//           Xs(j, col_counter) = X(i, j, s-1)*X(k, j, s-1);
//         }
//         col_counter++;
//       }
//     }
//   }
//
//
//   if(order == 3){
//     // third part (no process variables)
//     for(int i = 0; i < q-1; i++){
//       for(int k = i+1; k < q-2; k++){
//         for(int l = k+1; l < q; l++){
//           for(int j = 0; j < J; j++){
//             Xs(j, col_counter) = X(i, j, s-1)*X(k, j, s-1)*X(l, j, s-1);
//           }
//           col_counter++;
//         }
//       }
//     }
//   }
//
//
//
//   if(order == 4){
//     // Third part of the sum (with process variables)
//     for(int i = 0; i < n_pv; i++){
//       for(int k = 0; k < q; k++){
//
//         for(int j = 0; j < J; j++){
//           Xs(j, col_counter) = X(k, j, s-1)*X(i+q, j, s-1); //x_k*z_i
//         }
//         col_counter++;
//       }
//     }
//
//     // Fourth part of the sum (with process variables)
//     for(int i = 0; i < n_pv-1; i++){
//       for(int k = i+1; k < n_pv; k++){
//
//         for(int j = 0; j < J; j++){
//           Xs(j, col_counter) = X(i+q, j, s-1)*X(k+q, j, s-1); //z_i*z_k
//         }
//         col_counter++;
//       }
//     }
//
//     // Fifth part of the sum (with process variables)
//     for(int i = 0; i < n_pv; i++){
//
//       for(int j = 0; j < J; j++){
//         Xs(j, col_counter) = X(i+q, j, s-1)*X(i+q, j, s-1); //z_i^2
//       }
//       col_counter++;
//     }
//   }
//
//   return Xs;
// }




// [[Rcpp::export]]
arma::vec getUsMNL(arma::cube& X, arma::vec& beta, int s, arma::mat& Xs, bool transform_beta = true, int n_pv = 0){
  // Function that returns the utility vector of choice set s.
  // Final vector is of length J.
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
  // Input:
  //     X: design cube of dimensions (q, J, S)
  //     beta: parameter vector. Must be of length m or length m-1
  //     s: integer s, corresponding to a choice set in 1 to S.
  //     Xs: design matrix of choice set s. Must be of dimension (J, m-1).
  //     transform_beta: boolean parameter. Should the beta vector/matrix be transformed by subtracting the q-th element?

  int J = X.n_cols;
  int q = X.n_rows - n_pv;

  int m = Xs.n_cols + 1;

  // Check input dimensions
  // if(m != beta.n_elem) stop("Incompatible q in beta and X");

  arma::vec Us(J);

  if(transform_beta){
    // Create auxiliary vector
    arma::vec beta2(m-1);

    // compute beta_i_star = beta_i - beta_q
    for(int i = 0; i < q-1; i++){
      beta2(i) = beta(i) - beta(q-1);
    }

    for(int i = q-1; i < m-1; i++){
      beta2(i) = beta(i+1);
    }
    Us = Xs*beta2;
  } else{
    Us = Xs*beta;
  }

  return Us;
}


// [[Rcpp::export]]
arma::vec getPsMNL(arma::cube& X, arma::vec& beta, int s, arma::mat& Xs, bool transform_beta = true, int n_pv = 0){
  // Function that returns the probability vector of choice set s, based on the softmax function.
  // Final vector is of length J.

  int J = X.n_cols;

  arma::vec Us(J);
  Us = getUsMNL(X, beta, s, Xs, transform_beta, n_pv);

  arma::vec exp_Ujs(J);
  arma::vec P(J);

  // subtracting the maximum value to avoid numerical overflow
  exp_Ujs = exp(Us - max(Us));

  double sum_exp_Ujs = sum(exp_Ujs);

  P = exp_Ujs/sum_exp_Ujs;

  // Check for numerical instabilities
  if(abs(sum(P) - 1) > 1e-10) warning("Sum may not be numerically equal to 1.");

  return P;
}





// [[Rcpp::export]]
arma::mat getInformationMatrixMNL(arma::cube& X, arma::vec& beta, int order, bool transform_beta = true, int n_pv = 0, bool no_choice = false){
  // Function that returns the information matrix for design cube X and parameter vector beta.
  // It is the sum of the information matrices of the S choice sets.
  // For more information about this function, see R's wrapper function mnl_get_information_matrix().

  int J = X.n_cols;
  int S = X.n_elem/(X.n_cols*X.n_rows);

  // Get m
  // Should probably write this in a function
  int q = X.n_rows - n_pv;
  int m;
  if(order == 1){
    m = q;
  } else{
    if(order == 2){
      m = q*(q-1)/2 + q;
    } else{
      if(order == 3){
        m = (q*q*q + 5*q)/6; // = q + q*(q-1)/2 + q*(q-1)*(q-2)/6 = (q^3+ 5*q)/6
      } else{
        if(order == 4){
          m = q + (q-1)*q/2 + q*n_pv + n_pv*(n_pv-1)/2 + n_pv;
        } else{
          // Error. Not necessary to check here because there are input checks in R.
          stop("Wrong order.");
        }
      }
    }
  }

  // if(!transform_beta) m = m-1;

  // int m;
  // if(transform_beta) m = beta.n_elem;
  // else m = beta.n_elem + 1;


  int nrow_Xs, ncol_Xs;
  if(!no_choice){
    nrow_Xs = J;
    ncol_Xs = m-1;
  } else{
    nrow_Xs = J+1;
    ncol_Xs = m;
  }

  arma::mat Xs(nrow_Xs, ncol_Xs);

  arma::mat I(ncol_Xs, ncol_Xs, fill::zeros);
  // arma::mat identity(J, J, fill::eye);
  arma::mat ps_ps_t(J, J);
  arma::mat middle(J, J);
  arma::vec ps;

  // Compute information matrix for each choice set s, and sum.
  for(int s = 1; s <= S; s++){
    Xs = getXsMNL(X, s, order, n_pv, no_choice);
    ps = getPsMNL(X, beta, s, Xs, transform_beta, n_pv);
    ps_ps_t = ps*ps.t();
    middle = ps_ps_t;
    middle.diag() = ps_ps_t.diag() - ps;
    I = I - (Xs.t())*middle*Xs;
  }
  return I;
}


arma::mat getInformationMatrixMNL(arma::cube& X, arma::cube& Xs_cube, arma::vec& beta, int order, bool transform_beta = true, int n_pv = 0, bool no_choice = false){
  // Redefinition of getInformationMatrixMNL() so that it has an extra parameter, Xs_cube. This extra parameter is an Armadillo cube that
  // has the Xs matrices precomputed for each choice set s. It is useful to save time in the Bayesian case so that the matrix Xs is not
  // recomputed with each draw from the prior distribution.

  int J = X.n_cols;
  int S = X.n_elem/(X.n_cols*X.n_rows);

  // Get m
  // Should probably write this in a function
  int q = X.n_rows - n_pv;
  int m;
  if(order == 1){
    m = q;
  } else{
    if(order == 2){
      m = q*(q-1)/2 + q;
    } else{
      if(order == 3){
        m = (q*q*q + 5*q)/6; // = q + q*(q-1)/2 + q*(q-1)*(q-2)/6 = (q^3+ 5*q)/6
      } else{
        if(order == 4){
          m = q + (q-1)*q/2 + q*n_pv + n_pv*(n_pv-1)/2 + n_pv;
        } else{
          // Error. Not necessary to check here because there are input checks in R.
          stop("Wrong order.");
        }
      }
    }
  }

  // if(!transform_beta) m = m-1;

  // int m;
  // if(transform_beta) m = beta.n_elem;
  // else m = beta.n_elem + 1;


  int nrow_Xs, ncol_Xs;
  if(!no_choice){
    nrow_Xs = J;
    ncol_Xs = m-1;
  } else{
    nrow_Xs = J+1;
    ncol_Xs = m;
  }

  arma::mat Xs(nrow_Xs, ncol_Xs);

  arma::mat I(ncol_Xs, ncol_Xs, fill::zeros);
  // arma::mat identity(J, J, fill::eye);
  arma::mat ps_ps_t(J, J);
  arma::mat middle(J, J);
  arma::vec ps;

  // Compute information matrix for each choice set s, and sum.
  for(int s = 1; s <= S; s++){
    Xs = Xs_cube(arma::span(s-1), arma::span::all, arma::span::all);
    ps = getPsMNL(X, beta, s, Xs, transform_beta, n_pv);
    ps_ps_t = ps*ps.t();
    middle = ps_ps_t;
    middle.diag() = ps_ps_t.diag() - ps;
    I = I - (Xs.t())*middle*Xs;
  }
  return I;
}






// [[Rcpp::export]]
arma::mat getBayesianInformationMatrixMNL(arma::cube& X, arma::mat& beta_mat, int order, bool transform_beta = true, int n_pv = 0, bool no_choice = false){
  // Computes optimality value.
  // For more information see R's wrapper function mnl_get_opt_crit_value().


  int J = X.n_cols;
  int S = X.n_elem/(X.n_cols*X.n_rows);

  // Get m
  // Should probably write this in a function
  int q = X.n_rows - n_pv;
  int m;
  if(order == 1){
    m = q;
  } else{
    if(order == 2){
      m = q*(q-1)/2 + q;
    } else{
      if(order == 3){
        m = (q*q*q + 5*q)/6; // = q + q*(q-1)/2 + q*(q-1)*(q-2)/6 = (q^3+ 5*q)/6
      } else{
        if(order == 4){
          m = q + (q-1)*q/2 + q*n_pv + n_pv*(n_pv-1)/2 + n_pv;
        } else{
          // Error. Not necessary to check here because there are input checks in R.
          stop("Wrong order.");
        }
      }
    }
  }

  // if(!transform_beta) m = m-1;

  // int m;
  // if(transform_beta) m = beta.n_elem;
  // else m = beta.n_elem + 1;


  int nrow_Xs, ncol_Xs;
  if(!no_choice){
    nrow_Xs = J;
    ncol_Xs = m-1;
  } else{
    nrow_Xs = J+1;
    ncol_Xs = m;
  }

  arma::mat Xs(nrow_Xs, ncol_Xs);

  arma::mat I(ncol_Xs, ncol_Xs, fill::zeros);


  // int m;
  // if(transform_beta) m = beta_mat.n_cols;
  // else m = beta_mat.n_cols + 1;
  //
  // // Initialize information matrix with zeros
  // arma::mat I(m-1, m-1, fill::zeros);
  // vec beta = zeros(m-1);

  vec beta = zeros(ncol_Xs);
  int n_sims = beta_mat.n_rows;

  // Create the cube Xs_cube with the Xs matrix for each choice set s.
  // arma::mat Xs(J, ncol_Xs);
  arma::cube Xs_cube(S, nrow_Xs, ncol_Xs);
  for(int s = 1; s <= S; s++){
    Xs = getXsMNL(X, s, order, n_pv, no_choice);
    // Copy the information from Xs to the cube. Still can't find a way to copy a matrix to a subcube in Armadillo.
    for(int ix1 = 0; ix1 < nrow_Xs; ix1++){
      for(int ix2 = 0; ix2 < ncol_Xs; ix2++){
        Xs_cube(s-1, ix1, ix2) = Xs(ix1, ix2);
      }
    }
  }

  // Iterate over all prior draws
  for(int i = 0; i < n_sims; i++){
    beta = conv_to<vec>::from(beta_mat.row(i));
    I += getInformationMatrixMNL(X, Xs_cube, beta, order, transform_beta, n_pv, no_choice);

  } // end for


  return I/n_sims;
}




// [[Rcpp::export]]
double getOptCritValueMNL(arma::cube& X, arma::mat& beta_mat, int verbose, int opt_crit, arma::mat& W, int order, bool transform_beta = true, int n_pv = 0, bool no_choice = false){
  // Computes optimality value.
  // For more information see R's wrapper function mnl_get_opt_crit_value().


  int J = X.n_cols;
  int S = X.n_elem/(X.n_cols*X.n_rows);

  // Get m
  // Should probably write this in a function
  int q = X.n_rows - n_pv;
  int m;
  if(order == 1){
    m = q;
  } else{
    if(order == 2){
      m = q*(q-1)/2 + q;
    } else{
      if(order == 3){
        m = (q*q*q + 5*q)/6; // = q + q*(q-1)/2 + q*(q-1)*(q-2)/6 = (q^3+ 5*q)/6
      } else{
        if(order == 4){
          m = q + (q-1)*q/2 + q*n_pv + n_pv*(n_pv-1)/2 + n_pv;
        } else{
          // Error. Not necessary to check here because there are input checks in R.
          stop("Wrong order.");
        }
      }
    }
  }

  // if(!transform_beta) m = m-1;

  // int m;
  // if(transform_beta) m = beta.n_elem;
  // else m = beta.n_elem + 1;


  int nrow_Xs, ncol_Xs;
  if(!no_choice){
    nrow_Xs = J;
    ncol_Xs = m-1;
  } else{
    nrow_Xs = J+1;
    ncol_Xs = m;
  }

  arma::mat Xs(nrow_Xs, ncol_Xs);

  arma::mat I(ncol_Xs, ncol_Xs, fill::zeros);


  // int m;
  // if(transform_beta) m = beta_mat.n_cols;
  // else m = beta_mat.n_cols + 1;
  //
  // // Initialize information matrix with zeros
  // arma::mat I(m-1, m-1, fill::zeros);
  // vec beta = zeros(m-1);

  vec beta = zeros(ncol_Xs);
  int n_sims = beta_mat.n_rows;

  // Final efficiency criterion value. We want to minimize this.
  double eff_crit_val = 0.0;

  // Accumulator
  double acc = 0.0;



  bool opt_success; // flag for success of the decomposition



  // Create the cube Xs_cube with the Xs matrix for each choice set s.
  // arma::mat Xs(J, ncol_Xs);
  arma::cube Xs_cube(S, nrow_Xs, ncol_Xs);
  for(int s = 1; s <= S; s++){
    Xs = getXsMNL(X, s, order, n_pv, no_choice);
    // Copy the information from Xs to the cube. Still can't find a way to copy a matrix to a subcube in Armadillo.
    for(int ix1 = 0; ix1 < nrow_Xs; ix1++){
      for(int ix2 = 0; ix2 < ncol_Xs; ix2++){
        Xs_cube(s-1, ix1, ix2) = Xs(ix1, ix2);
      }
    }
  }

  // Iterate over all prior draws
  for(int i = 0; i < n_sims; i++){
    beta = conv_to<vec>::from(beta_mat.row(i));
    I = getInformationMatrixMNL(X, Xs_cube, beta, order, transform_beta, n_pv, no_choice);
    if(verbose >= 5) Rcout << "Information matrix. I = \n" << I << std::endl;

    if(opt_crit == 0){
      // D-optimality
      arma::mat L; // matrix for Cholesky decomposition

      opt_success = chol(L, I);
      // If the decomposition fails chol(R,X) resets R and returns a bool set to false (exception is not thrown) (http://arma.sourceforge.net/docs.html#chol)
      if(opt_success){
        acc = acc - 2*sum(log(L.diag()));
      } else {
        // If Cholesky decomposition fails, it is likely because information matrix was not numerically positive definite.
        break;
      }

    } else{
      // I-optimality
      arma::mat C; // Matrix for the solution of the linear system involving I and W.

      opt_success = solve(C, I, W, solve_opts::likely_sympd + solve_opts::no_approx);
      // If no solution is found solve(X,A,B) resets X and returns a bool set to false (exception is not thrown)
      if(opt_success){
        acc = acc + log(trace(C));
      } else {
        break;
      }
    }

  } // end for

  if(!opt_success){
    // If either the Cholesky decomposition or the solving of the linear system fail, it is probably because a numerical inestability.
    // The function then returns the efficiency value as a big positive number (1000), this way the algorithm does nothing in this iteration
    // because it thinks there was no improvement when swapping the proportions.
    eff_crit_val = 1000;
    warning("Numerical problem in Cholesky decomposition or while solving linear system.");
  } else{
    eff_crit_val = acc/n_sims;
  }

  return eff_crit_val;
}







// [[Rcpp::export]]
void findBestCoxDirMNLDiscrete( // void but modifies X
    arma::mat& cox_dir, arma::cube& X, arma::mat& beta_mat,
    int k, int s, double opt_crit_value_best,
    int verbose, int opt_crit, arma::mat& W, int order, bool transform_beta = true) {
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

    opt_crit_value_j = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W, order, transform_beta);

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
void changeIngredientDesignMNL(double theta, arma::cube& X, int i, int j, int s, int n_pv){ // Void but modifies X
  // Modifies cube X changing the i-th ingredient in the j-th alternative in s-th choice set to theta.
  // Since theta is an ingredient proportion, it must be between 0 and 1.
  // Indices i, j and k are 0-indexed.


  // Number of ingredients
  int q = X.n_rows - n_pv;

  // Create a vector with the ingredients of the j-th alternative of the s-th choice set.
  arma::vec x = X(arma::span(0, q-1), arma::span(j), arma::span(s));

  double delta = theta - x(i);

  // recompute proportions:
  vec setDiff_aux = linspace<vec>(0, q-1, q);
  vec setDiff = removeElement(setDiff_aux, i);
  int k;

  for(int k_aux = 0; k_aux < setDiff.n_elem; k_aux++){
    k = setDiff(k_aux);

    if(abs(1 - x(i)) < 1e-13) { // In case x(i) is numerically 1
      // Two cases:
      // 1) If delta is numerically -1, it means that the change is from 1 to 0. Then, the rest of the ingredients must be 1/(q-1)
      // 2) If delta is not -1 and x(i) is numerically 1, then it means that the rest of the ingredients were 0, so they were in equal
      //    proportions and should remain that way.
      x(k) = (1.0 - theta)/((double)q - 1.0); // Same as x(k) = (-delta)/(q - 1.0); and x(k) = (1-delta-x(i))/(q - 1.0);
    } else{ // In case x(i) is not numerically 1
      // Other case
      x(k) = x(k) - delta*x(k)/(1 - x(i));
    }
  }



  x(i) = theta;

  if(abs(sum(x) - 1) > 1e-13){
    // Do not change design
    warning("Mixture ingredients do not sum up to numerical 1. Not changing this run of the design.");
    // Rcout << x << std::endl;
  } else{
    // Replace the design X with the recomputed proportions according to Cox direction
    for(int l = 0; l < q; l++){
      X(l, j, s) = x(l);
    }
  }



}




// [[Rcpp::export]]
double efficiencyCoxScheffeMNL(double theta, arma::cube& X, arma::mat& beta_mat,
                               int i, int j, int s,
                               int opt_crit, arma::mat& W, int order, bool transform_beta = true, int n_pv = 0, bool no_choice = false){
  // Computes efficiency criterion of a design cube X but where the i-th ingredient in the j-th
  // alternative in s-th choice set is changed to theta.
  // Since theta is an ingredient proportion, it must be between 0 and 1.
  // Indices i, j and k are 0-indexed.
  // We want to minimize this.

  int q = X.n_rows - n_pv;

  // Temporarily store the corresponding row in a vector
  arma::vec x_row(q);
  for(int l = 0; l < q; l++){
    x_row(l) = X(l, j, s);
  }

  changeIngredientDesignMNL(theta, X, i, j,s, n_pv);

  // Utility function value. We want to minimize this.
  double utility_funct_value = getOptCritValueMNL(X, beta_mat, 0, opt_crit, W, order, transform_beta, n_pv, no_choice);

  // return X to its original value
  for(int l = 0; l < q; l++){
    X(l, j, s) = x_row(l);
  }

  return(utility_funct_value);

}




double efficiencyPVScheffeMNL(double theta, arma::cube& X, arma::mat& beta_mat,
                              int i, int j, int s,
                              int opt_crit, arma::mat& W, int order, bool transform_beta = true, int n_pv = 0, bool no_choice = false){
  // Computes efficiency criterion of a design cube X but where the i-th process variable in the j-th
  // alternative in the s-th choice set is changed to theta.
  // Note: it should be shifted by the q proportion variables. That is, this function changes the i-th element
  //       without checking whether it is a process variable or an ingredient proportion. This is done in the main algorithm.
  // Here theta is NOT an ingredient proportion.
  // Indices i, j, and k are 0-indexed.
  // We want to minimize this.

  double original = X(i, j, s);

  // change value
  X(i, j, s) = theta;

  // Utility function value. We want to minimize this.
  double utility_funct_value = getOptCritValueMNL(X, beta_mat, 0, opt_crit, W, order, transform_beta, n_pv, no_choice);

  // return X to its original value
  X(i, j, s) = original;

  return(utility_funct_value);
}




void findBestPVMNLBrent( // void but modifies X
    arma::cube& X, arma::mat& beta_mat, int i, int j, int s, int opt_crit, int order, arma::mat& W,
    double lower = 0, double upper = 1, double tol = 0.0001, int verbose = 0, bool transform_beta = true, int n_pv = 0, bool no_choice = false) {
  // Must change upper and lower values

  // The optimality criterion value in the design that is being used as input
  double f_original = getOptCritValueMNL(X, beta_mat, 0, opt_crit, W, order, transform_beta, n_pv, no_choice);

  // Helper function that depends only on theta (the ingredient proportion that is being changed now)
  auto f = [&X, &beta_mat, i, j, s, opt_crit, &W, order, transform_beta, n_pv, no_choice](double theta){
    return efficiencyPVScheffeMNL(theta, X, beta_mat, i, j, s, opt_crit, W, order, transform_beta, n_pv, no_choice);
  };

  double theta_brent, theta_star, f_star;
  double f_brent = brent::local_min_mb(lower, upper, tol, f, theta_brent);

  // Check end points
  double f_lower = efficiencyPVScheffeMNL(lower, X, beta_mat, i, j, s, opt_crit, W, order, transform_beta, n_pv, no_choice);
  double f_upper = efficiencyPVScheffeMNL(upper, X, beta_mat, i, j, s, opt_crit, W, order, transform_beta, n_pv, no_choice);

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
    X(i, j, s) = theta_star;
  }
}




// [[Rcpp::export]]
void findBestCoxDirMNLBrent( // void but modifies X
    arma::cube& X, arma::mat& beta_mat, int i, int j, int s, int opt_crit, int order, arma::mat& W,
    double lower = -1, double upper = 1, double tol = 0.0001, int verbose = 0, bool transform_beta = true, int n_pv = 0, bool no_choice = false) {
  // Finds the best point in the Cox direction using Brent's optimization method


  // The optimality criterion value in the design that is being used as input
  double f_original = getOptCritValueMNL(X, beta_mat, 0, opt_crit, W, order, transform_beta, n_pv, no_choice);


  // Helper function that depends only on theta (the ingredient proportion that is being changed now)
  auto f = [&X, &beta_mat, i, j, s, opt_crit, &W, order, transform_beta, n_pv, no_choice](double theta){
    return efficiencyCoxScheffeMNL(theta, X, beta_mat, i, j, s, opt_crit, W, order, transform_beta, n_pv, no_choice);
  };

  // theta_brent: the "optimal" theta that is going to be returned by Brent's method
  // theta_star: the actual optimal theta
  // f_star: the value of the objective function evaluated in the optimal theta
  double theta_brent, theta_star, f_star;
  double f_brent = brent::local_min_mb(lower, upper, tol, f, theta_brent);

  // Check end points:
  double f_lower = efficiencyCoxScheffeMNL(lower, X, beta_mat, i, j, s, opt_crit, W, order, transform_beta, n_pv, no_choice);
  double f_upper = efficiencyCoxScheffeMNL(upper, X, beta_mat, i, j, s, opt_crit, W, order, transform_beta, n_pv, no_choice);

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

  if(verbose >= 5){
    Rcout << "\ttheta_star: " << theta_star << "\tf_star: " << f_star << std::endl;
  }

  // If Brent's method didn't do any improvement, leave the original design as it is.
  // If it did, replace the design with theta star
  if(f_original > f_star){
    changeIngredientDesignMNL(theta_star, X, i, j, s, n_pv);
  }
}












// [[Rcpp::export]]
Rcpp::List mixtureCoordinateExchangeMNL(
    arma::cube X_orig, arma::mat beta_mat, int order, int max_it, int verbose, int opt_crit, arma::mat W,
    int opt_method, double lower, double upper, double tol, int n_cox_points, bool transform_beta = true, int n_pv = 0, bool no_choice = false){
  // See mnl_mixture_coord_exch() in R for details.
  // lower: lower bound in Brent's optimization method
  // upper: upper bound in Brent's optimization method
  // W: moment matrix


  // Does not do input checks because the R wrapper function does them.
  // (Except for opt method. Have to remove later and put it in R.)
  if(opt_method == 1 & n_pv > 0) {
    stop("Discretization of Cox direction is not available for models that include process variables.");
  }

  // Create a vector to store the values of the efficiency metric in each iteration.
  arma::vec efficiency_value_per_iteration(max_it + 1, fill::zeros);

  // Create new cube, otherwise it is modified in R too
  arma::cube X = X_orig;

  int J = X.n_cols;
  int q = X.n_rows - n_pv;
  int S = X.n_elem/(X.n_cols*X.n_rows);

  // Create matrix with appropriate dimensions for Cox direction in each iteration
  arma::mat cox_dir(n_cox_points, q);

  // Vector of ingredient proportions
  arma::vec x(q);

  double opt_crit_value_orig = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W, order, transform_beta, n_pv, no_choice);
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
            findBestCoxDirMNLBrent(X, beta_mat, i, k-1, s-1, opt_crit, order, W, lower, upper, tol, verbose, transform_beta, n_pv, no_choice);
          } else{
            if(no_choice) stop("Cannot use discretization when no_choice = T");
            cox_dir = computeCoxDirection(x, i+1, n_cox_points, verbose);
            findBestCoxDirMNLDiscrete(cox_dir, X, beta_mat, k, s, opt_crit_value_best, verbose, opt_crit, W, order, transform_beta);
          }

          opt_crit_value_best = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W, order, transform_beta, n_pv, no_choice);

          if(verbose >= 2) Rcout << "Opt-crit-value: " << opt_crit_value_best << std::endl;

          if(verbose >= 5){
            Rcout << "X =\n" << X << std::endl;
          }

        } // end for i



        // Process variables
        if(n_pv > 0){
          for(int pv = q; pv < n_pv + q; pv++){

            if(verbose >= 2) Rcout << "\nIter: " << it <<  ", k = " << k << ", s = " << s << ", pv = " << pv << std::endl;

            // findBestCoxDirMNLBrent(X, beta_mat, i, k-1, s-1, opt_crit, order, W, lower, upper, tol, verbose, transform_beta, n_pv);
            // (X, run-1, pv, order, opt_crit, W, lower, upper, tol, n_pv);
            findBestPVMNLBrent(X, beta_mat, pv, k-1, s-1, opt_crit, order, W, -1, 1, tol, verbose, transform_beta, n_pv, no_choice);

            opt_crit_value_best = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W, order, transform_beta, n_pv, no_choice);

            if(verbose >= 2) Rcout << "Opt-crit-value: " << opt_crit_value_best << std::endl;

            if(verbose >= 5){
              Rcout << "X =\n" << X << std::endl;
            }

          }
        }



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





















// [[Rcpp::export]]
void changeIngredientDesignMNL2(double theta, arma::cube& X, int i, arma::vec j_vec, arma::vec s_vec, int n_pv){ // Void but modifies X
  // Modifies cube X changing the i-th ingredient in alternatives j_vec and choice sets s_vec to theta.
  // Since theta is an ingredient proportion, it must be between 0 and 1.
  // Indices i, j_vec and s_vec are 0-indexed.

  int j_len = j_vec.n_elem;
  int s_len = s_vec.n_elem;

  if(j_len != s_len) stop("j_len and s_len should have same length");

  int j;
  int s;


  // Number of ingredients
  int q = X.n_rows - n_pv;

  // Create a vector with the ingredients of the j-th alternative of the s-th choice set.
  arma::vec x = X(arma::span(0, q-1), arma::span(j_vec(0)), arma::span(s_vec(0)));

  double delta = theta - x(i);

  // recompute proportions:
  vec setDiff_aux = linspace<vec>(0, q-1, q);
  vec setDiff = removeElement(setDiff_aux, i);
  int k;

  for(int k_aux = 0; k_aux < setDiff.n_elem; k_aux++){
    k = setDiff(k_aux);

    if(abs(1 - x(i)) < 1e-13) { // In case x(i) is numerically 1
      // Two cases:
      // 1) If delta is numerically -1, it means that the change is from 1 to 0. Then, the rest of the ingredients must be 1/(q-1)
      // 2) If delta is not -1 and x(i) is numerically 1, then it means that the rest of the ingredients were 0, so they were in equal
      //    proportions and should remain that way.
      x(k) = (1.0 - theta)/((double)q - 1.0); // Same as x(k) = (-delta)/(q - 1.0); and x(k) = (1-delta-x(i))/(q - 1.0);
    } else{ // In case x(i) is not numerically 1
      // Other case
      x(k) = x(k) - delta*x(k)/(1 - x(i));
    }
  }


  x(i) = theta;

  if(abs(sum(x) - 1) > 1e-13){
    // Do not change design
    warning("Mixture ingredients do not sum up to numerical 1. Not changing this run of the design.");
    Rcout << x << std::endl;
  } else{
    // Replace the design X with the recomputed proportions according to Cox direction
    for(int l = 0; l < q; l++){
      for(int i1 = 0; i1 < j_len; i1++){
        j = j_vec(i1);
        s = s_vec(i1);
        X(l, j, s) = x(l);
      }

    }
  }
}




// [[Rcpp::export]]
double efficiencyCoxScheffeMNL2(double theta, arma::cube& X, arma::mat& beta_mat,
                                int i, arma::vec j_vec, arma::vec s_vec,
                                int opt_crit, arma::mat& W, int order, bool transform_beta = true, int n_pv = 0, bool no_choice = false){
  // Computes efficiency criterion of a design cube X but where the i-th ingredient in the j-th
  // alternative in s-th choice set is changed to theta.
  // Since theta is an ingredient proportion, it must be between 0 and 1.
  // Index vectors j_vec and s_vec, as well as integer i, are 0-indexed.
  // We want to minimize this.


  int j_len = j_vec.n_elem;
  int s_len = s_vec.n_elem;

  if(j_len != s_len) stop("j_len and s_len should have same length");

  int j;
  int s;


  int q = X.n_rows - n_pv;

  // Temporarily store the corresponding row in a vector
  arma::vec x_row(q);
  for(int l = 0; l < q; l++){
    x_row(l) = X(l, j_vec(0), s_vec(0));
  }

  changeIngredientDesignMNL2(theta, X, i, j_vec, s_vec, n_pv);

  // Utility function value. We want to minimize this.
  double utility_funct_value = getOptCritValueMNL(X, beta_mat, 0, opt_crit, W, order, transform_beta, n_pv, no_choice);

  // return X to its original value
  for(int l = 0; l < q; l++){
    for(int i1 = 0; i1 < j_len; i1++){
      j = j_vec(i1);
      s = s_vec(i1);
      X(l, j, s) = x_row(l);
    }
  }

  return(utility_funct_value);

}








void findBestCoxDirMNLBrent2( // void but modifies X
    arma::cube& X, arma::mat& beta_mat, int i, arma::vec j_vec, arma::vec s_vec, int opt_crit, int order, arma::mat& W,
    double lower = -1, double upper = 1, double tol = 0.0001, int verbose = 0, bool transform_beta = true, int n_pv = 0, bool no_choice = false) {
  // Finds the best point in the Cox direction using Brent's optimization method


  // The optimality criterion value in the design that is being used as input
  double f_original = getOptCritValueMNL(X, beta_mat, 0, opt_crit, W, order, transform_beta, n_pv, no_choice);


  // Helper function that depends only on theta (the ingredient proportion that is being changed now)
  auto f = [&X, &beta_mat, i, j_vec, s_vec, opt_crit, &W, order, transform_beta, n_pv, no_choice](double theta){
    return efficiencyCoxScheffeMNL2(theta, X, beta_mat, i, j_vec, s_vec, opt_crit, W, order, transform_beta, n_pv, no_choice);
  };

  // theta_brent: the "optimal" theta that is going to be returned by Brent's method
  // theta_star: the actual optimal theta
  // f_star: the value of the objective function evaluated in the optimal theta
  double theta_brent, theta_star, f_star;
  double f_brent = brent::local_min_mb(lower, upper, tol, f, theta_brent);

  // Check end points:
  double f_lower = efficiencyCoxScheffeMNL2(lower, X, beta_mat, i, j_vec, s_vec, opt_crit, W, order, transform_beta, n_pv, no_choice);
  double f_upper = efficiencyCoxScheffeMNL2(upper, X, beta_mat, i, j_vec, s_vec, opt_crit, W, order, transform_beta, n_pv, no_choice);

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

  if(verbose >= 5){
    Rcout << "\ttheta_star: " << theta_star << "\tf_star: " << f_star << std::endl;
  }

  // If Brent's method didn't do any improvement, leave the original design as it is.
  // If it did, replace the design with theta star
  if(f_original > f_star){
    changeIngredientDesignMNL2(theta_star, X, i, j_vec, s_vec, n_pv);
  }
}












// Needs to be updated to use no_choice
// [[Rcpp::export]]
Rcpp::List mixtureCoordinateExchangeMNL2(
    arma::cube X_orig, arma::mat beta_mat, int order, int max_it, int verbose, int opt_crit, arma::mat W,
    int opt_method, double lower, double upper, double tol, int n_cox_points, bool transform_beta, int n_pv, bool no_choice,
    arma::mat mapping_matrix){
  // j_vec and s_vec are 1-indexed

  bool do_clustering = mapping_matrix.n_elem > 1;

  arma::vec s_vec_all;
  arma::vec j_vec_all;
  arma::vec point_id_vec;
  arma::vec s_vec;
  arma::vec j_vec;
  arma::vec s_vec_one_element(1);
  arma::vec j_vec_one_element(1);
  if(do_clustering){
    s_vec_all = mapping_matrix(span::all, 0);
    j_vec_all = mapping_matrix(span::all, 1);
    point_id_vec = mapping_matrix(span::all, 2);
  }


  // See mnl_mixture_coord_exch() in R for details.
  // lower: lower bound in Brent's optimization method
  // upper: upper bound in Brent's optimization method
  // W: moment matrix


  // Does not do input checks because the R wrapper function does them.
  // (Except for opt method. Have to remove later and put it in R.)
  if(opt_method == 1 & n_pv > 0) {
    stop("Discretization of Cox direction is not available for models that include process variables.");
  }

  // Create a vector to store the values of the efficiency metric in each iteration.
  arma::vec efficiency_value_per_iteration(max_it + 1, fill::zeros);

  // Create new cube, otherwise it is modified in R too
  arma::cube X = X_orig;

  int J = X.n_cols;
  int q = X.n_rows - n_pv;
  int S = X.n_elem/(X.n_cols*X.n_rows);

  // Create matrix with appropriate dimensions for Cox direction in each iteration
  arma::mat cox_dir(n_cox_points, q);

  // Vector of ingredient proportions
  arma::vec x(q);

  double opt_crit_value_orig = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W, order, transform_beta, n_pv, no_choice);
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

    int row_counter = 0;

    for(int j = 1; j <= J; j++){

      for(int s = 1; s <= S; s++){

        if(do_clustering){
          int point_id = mapping_matrix(row_counter, 2);
          uvec point_id_ix = find(point_id_vec == point_id);
          s_vec = s_vec_all(point_id_ix) - 1;
          j_vec = j_vec_all(point_id_ix) - 1;
        } else{
          s_vec_one_element(0) = s - 1;
          j_vec_one_element(0) = j - 1;

        }


        // Rcout << "row_counter: " << row_counter << std::endl;
        // Rcout << s_vec << std::endl;
        // Rcout << j_vec << std::endl;

        for(int i = 0; i < q; i++){
          // if(verbose >= 2) Rcout << "\t\ti = " << i << std::endl;
          if(verbose >= 2) Rcout << "\nIter: " << it <<  ", j = " << j << ", s = " << s << ", i = " << i << std::endl;

          // populate x vector with corresponding ingredient proportions
          for(int l = 0; l < q; l++){
            x(l) = X(l, j-1, s-1);
          }

          if(opt_method == 0){
            if(do_clustering){
              findBestCoxDirMNLBrent2(X, beta_mat, i, j_vec, s_vec, opt_crit, order, W, lower, upper, tol, verbose, transform_beta, n_pv, no_choice);
            } else{
              findBestCoxDirMNLBrent2(X, beta_mat, i, j_vec_one_element, s_vec_one_element, opt_crit, order, W, lower, upper, tol, verbose, transform_beta, n_pv, no_choice);
            }
          } else{
            cox_dir = computeCoxDirection(x, i+1, n_cox_points, verbose);
            findBestCoxDirMNLDiscrete(cox_dir, X, beta_mat, j, s, opt_crit_value_best, verbose, opt_crit, W, order, transform_beta);
          }

          opt_crit_value_best = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W, order, transform_beta, n_pv, no_choice);

          if(verbose >= 2) Rcout << "Opt-crit-value: " << opt_crit_value_best << std::endl;

          if(verbose >= 5){
            Rcout << "X =\n" << X << std::endl;
          }

        } // end for i



        // Process variables
        if(n_pv > 0){
          for(int pv = q; pv < n_pv + q; pv++){

            if(verbose >= 2) Rcout << "\nIter: " << it <<  ", j = " << j << ", s = " << s << ", pv = " << pv << std::endl;

            // findBestCoxDirMNLBrent(X, beta_mat, i, j-1, s-1, opt_crit, order, W, lower, upper, tol, verbose, transform_beta, n_pv);
            // (X, run-1, pv, order, opt_crit, W, lower, upper, tol, n_pv);
            findBestPVMNLBrent(X, beta_mat, pv, j-1, s-1, opt_crit, order, W, -1, 1, tol, verbose, transform_beta, n_pv, no_choice);

            opt_crit_value_best = getOptCritValueMNL(X, beta_mat, verbose, opt_crit, W, order, transform_beta, n_pv, no_choice);

            if(verbose >= 2) Rcout << "Opt-crit-value: " << opt_crit_value_best << std::endl;

            if(verbose >= 5){
              Rcout << "X =\n" << X << std::endl;
            }

          }
        }

        row_counter++;

      } // end for s

    } // end for j


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







