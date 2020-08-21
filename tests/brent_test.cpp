#include <Rcpp.h>
#include "brent.hpp"
using namespace Rcpp;

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
  double f_star = brent::local_min(-10.0, 10.0, 0.0001, my_f2, x_star);
  return List::create(Named("minimum") = x_star, Named("objective") = f_star);
}


// [[Rcpp::export]]
List min_f2_2(){
  double x_star;
  double f_star = brent::local_min(-10.0, 10.0, 0.0001, f2, x_star);
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
  double f_star = brent::local_min(-10.0, 10.0, 0.0001, f3, x_star);
  return List::create(Named("minimum") = x_star, Named("objective") = f_star);
}




// [[Rcpp::export]]
double banana_xy(double x, double y){
  return (1 - x)*(1-x) + 100*(y - x*x)*(y - x*x);
}


// [[Rcpp::export]]
double banana_0(double x, double y){
  return banana_xy(x, 1.0);
}


// [[Rcpp::export]]
double banana(double x){
  double y = 1.0;
  return (1 - x)*(1-x) + 100*(y - x*x)*(y - x*x);
}


// [[Rcpp::export]]
List min_banana(double lwr_bnd = -10, double uppr_bnd = 10, double tol = 0.0001){
  double x_star;
  double f_star = brent::local_min(lwr_bnd, uppr_bnd, tol, banana, x_star);
  return List::create(Named("minimizer") = x_star, Named("objective_func") = f_star);
}
















