#include <RcppArmadillo.h>
#include <functional>
#include "brent.hpp"
#include "utils.hpp"


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















