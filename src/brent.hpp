#include <vector>
#include <armadillo>


namespace brent {

  class func_base{
    public:
      virtual double operator() (double) = 0;
  };

  // double local_min ( double a, double b, double t, func_base& f,
  //   double &x );

  double local_min_mb ( double a, double b, double t, std::function<double (double)> f, double &x);

  double r8_epsilon ( );
  double r8_max ( double x, double y );
  double r8_sign ( double x );
  void timestamp ( );

  // // === simple wrapper functions
  // // === for convenience and/or compatibility
  // double local_min ( double a, double b, double t, double f ( double x ),
  //   double &x );

}

