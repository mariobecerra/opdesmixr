#include <Rcpp.h>
#include "brent.hpp"
#include <functional>

using namespace Rcpp;
using namespace std;

class myFunctorClass2 : public brent::func_base
{
  private:
    double a, b;
  public:
    myFunctorClass2 (double a_, double b_) : a(a_), b(b_) {}
  double operator() (double x) {
    return (x - a) * (x - a) + b;
  }
};

// [[Rcpp::export]]
double f2(double x){
  myFunctorClass2 my_f22(0.0, 0.0);
  return my_f22(x);
}

// [[Rcpp::export]]
List min_f2(){
  myFunctorClass2 my_f2(0.0, 0.0);
  double x_star;
  double f_star = brent::local_min_mb(-10.0, 10.0, 0.0001, my_f2, x_star);
  return List::create(Named("minimum") = x_star, Named("objective") = f_star);
}


// [[Rcpp::export]]
List min_f2_2(){
  double x_star;
  double f_star = brent::local_min_mb(-10.0, 10.0, 0.0001, f2, x_star);
  return List::create(Named("minimum") = x_star, Named("objective") = f_star);
}


// [[Rcpp::export]]
double f3(double x){
  double a = 2;
  double b = 3;
  return (x - a) * (x - a) + b;
}


// [[Rcpp::export]]
List min_f3(){
  double x_star;
  double f_star = brent::local_min_mb(-10.0, 10.0, 0.0001, f3, x_star);
  return List::create(Named("minimum") = x_star, Named("objective") = f_star);
}




// [[Rcpp::export]]
double banana_xy(double x, double y){
  return (1 - x)*(1-x) + 100*(y - x*x)*(y - x*x);
}


// https://stackoverflow.com/questions/34994475/rcpp-function-to-construct-a-function
// https://github.com/NormaVTH
// https://github.com/manuelalcantara52
// https://codereview.stackexchange.com/questions/103762/implementation-of-brents-algorithm-to-find-roots-of-a-polynomial

std::function<double(double)> create_banana(double y) {
  auto out = [y](double x) {
    return(banana_xy(x, y));
    };
  return out;
}


// [[Rcpp::export]]
double banana_x_y1(double x){
  std::function<double(double)> foo = create_banana(1.0);
  return(foo(x));
}

// [[Rcpp::export]]
double banana_x_y2(double x){
  std::function<double(double)> foo = create_banana(2.0);
  return(foo(x));
}



// [[Rcpp::export]]
List min_banana_x_y1(double lwr_bnd = -10, double uppr_bnd = 10, double tol = 0.0001){
  double x_star;
  double f_star = brent::local_min_mb(lwr_bnd, uppr_bnd, tol, banana_x_y1, x_star);
  return List::create(Named("minimizer") = x_star, Named("objective_func") = f_star);
}





// [[Rcpp::export]]
List minimize_banana_fixed_y(double y = 1.0, double lwr_bnd = -10, double uppr_bnd = 10, double tol = 0.0001){
  double x_star;
  auto f = [y](double x){ return banana_xy(x, y); };

  double f_star = brent::local_min_mb(lwr_bnd, uppr_bnd, tol, f, x_star);
  return List::create(Named("minimizer") = x_star, Named("objective_func") = f_star);
}





// auto f = [](double x){ return (x+1) * (x+2) * (x+3); };
//
// double find_banana_minimizer(double y = 1.0, double lwr_bnd = -10, double uppr_bnd = 10, double tol = 0.0001){
//
//
//
//   double x_star;
//   brent::local_min_mb(lwr_bnd, uppr_bnd, tol, create_banana(y), x_star);
//   return x_star;
// }
















