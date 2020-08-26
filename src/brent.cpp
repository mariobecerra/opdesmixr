# include <cmath>
# include <cstdlib>
# include <ctime>
# include <iostream>

using namespace std;


# include "brent.hpp"


namespace brent{


  //****************************************************************************80

  double local_min_mb ( double a, double b, double t, std::function<double (double)> f, double &x)
  // Modification of local_min by Mario Becerra. Takes an std::function as input.

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
  //
  //  Discussion:
  //
  //    The method used is a combination of golden section search and
  //    successive parabolic interpolation.  Convergence is never much slower
  //    than that for a Fibonacci search.  If F has a continuous second
  //    derivative which is positive at the minimum (which is not at A or
  //    B), then convergence is superlinear, and usually of the order of
  //    about 1.324....
  //
  //    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
  //    F is never evaluated at two points closer than TOL.
  //
  //    If F is a unimodal function and the computed values of F are always
  //    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
  //    LOCAL_MIN approximates the abscissa of the global minimum of F on the
  //    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
  //
  //    If F is not unimodal, then LOCAL_MIN may approximate a local, but
  //    perhaps non-global, minimum to the same accuracy.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    17 July 2011
  //
  //  Author:
  //
  //    Original FORTRAN77 version by Richard Brent.
  //    C++ version by John Burkardt.
  //
  //  Reference:
  //
  //    Richard Brent,
  //    Algorithms for Minimization Without Derivatives,
  //    Dover, 2002,
  //    ISBN: 0-486-41998-3,
  //    LC: QA402.5.B74.
  //
  //  Parameters:
  //
  //    Input, double A, B, the endpoints of the interval.
  //
  //    Input, double T, a positive absolute error tolerance.
  //
  //    Input, func_base& F, a user-supplied c++ functor whose
  //    local minimum is being sought.  The input and output
  //    of F() are of type double.
  //
  //    Output, double &X, the estimated value of an abscissa
  //    for which F attains a local minimum value in [A,B].
  //
  //    Output, double LOCAL_MIN, the value F(X).
  //

  {
    double c;
    double d;
    double e;
    double eps;
    double fu;
    double fv;
    double fw;
    double fx;
    double m;
    double p;
    double q;
    double r;
    double sa;
    double sb;
    double t2;
    double tol;
    double u;
    double v;
    double w;
    //
    //  C is the square of the inverse of the golden ratio.
    //
    c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );

    eps = sqrt ( r8_epsilon ( ) );

    sa = a;
    sb = b;
    x = sa + c * ( b - a );
    w = x;
    v = w;
    e = 0.0;
    fx = f ( x );
    fw = fx;
    fv = fw;

    for ( ; ; )
    {
      m = 0.5 * ( sa + sb ) ;
      tol = eps * fabs ( x ) + t;
      t2 = 2.0 * tol;
      //
      //  Check the stopping criterion.
      //
      if ( fabs ( x - m ) <= t2 - 0.5 * ( sb - sa ) )
      {
        break;
      }
      //
      //  Fit a parabola.
      //
      r = 0.0;
      q = r;
      p = q;

      if ( tol < fabs ( e ) )
      {
        r = ( x - w ) * ( fx - fv );
        q = ( x - v ) * ( fx - fw );
        p = ( x - v ) * q - ( x - w ) * r;
        q = 2.0 * ( q - r );
        if ( 0.0 < q )
        {
          p = - p;
        }
        q = fabs ( q );
        r = e;
        e = d;
      }

      if ( fabs ( p ) < fabs ( 0.5 * q * r ) &&
           q * ( sa - x ) < p &&
           p < q * ( sb - x ) )
      {
        //
        //  Take the parabolic interpolation step.
        //
        d = p / q;
        u = x + d;
        //
        //  F must not be evaluated too close to A or B.
        //
        if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
        {
          if ( x < m )
          {
            d = tol;
          }
          else
          {
            d = - tol;
          }
        }
      }
      //
      //  A golden-section step.
      //
      else
      {
        if ( x < m )
        {
          e = sb - x;
        }
        else
        {
          e = sa - x;
        }
        d = c * e;
      }
      //
      //  F must not be evaluated too close to X.
      //
      if ( tol <= fabs ( d ) )
      {
        u = x + d;
      }
      else if ( 0.0 < d )
      {
        u = x + tol;
      }
      else
      {
        u = x - tol;
      }

      fu = f ( u );
      //
      //  Update A, B, V, W, and X.
      //
      if ( fu <= fx )
      {
        if ( u < x )
        {
          sb = x;
        }
        else
        {
          sa = x;
        }
        v = w;
        fv = fw;
        w = x;
        fw = fx;
        x = u;
        fx = fu;
      }
      else
      {
        if ( u < x )
        {
          sa = u;
        }
        else
        {
          sb = u;
        }

        if ( fu <= fw || w == x )
        {
          v = w;
          fv = fw;
          w = u;
          fw = fu;
        }
        else if ( fu <= fv || v == x || v == w )
        {
          v = u;
          fv = fu;
        }
      }
    }
    return fx;
  }




  // ======================================================================
  // === Simple wrapper functions
  // === for convenience and/or compatibility.
  //
  // === The three functions are the same as above,
  // === except that they take a plain function F
  // === instead of a c++ functor.  In all cases, the
  // === input and output of F() are of type double.

  typedef double DoubleOfDouble (double);

  class func_wrapper : public func_base {
    DoubleOfDouble* func;
  public:
    func_wrapper(DoubleOfDouble* f) {
      func = f;
    }
    virtual double operator() (double x){
      return func(x);
    }
  };




  double r8_epsilon ( )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    R8_EPSILON returns the R8 roundoff unit.
    //
    //  Discussion:
    //
    //    The roundoff unit is a number R which is a power of 2 with the
    //    property that, to the precision of the computer's arithmetic,
    //      1 < 1 + R
    //    but
    //      1 = ( 1 + R / 2 )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    01 September 2012
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Output, double R8_EPSILON, the R8 round-off unit.
    //
  {
    const double value = 2.220446049250313E-016;

    return value;
  }
  //****************************************************************************80

  double r8_max ( double x, double y )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    R8_MAX returns the maximum of two R8's.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    18 August 2004
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Input, double X, Y, the quantities to compare.
    //
    //    Output, double R8_MAX, the maximum of X and Y.
    //
  {
    double value;

    if ( y < x )
    {
      value = x;
    }
    else
    {
      value = y;
    }
    return value;
  }
  //****************************************************************************80

  double r8_sign ( double x )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    R8_SIGN returns the sign of an R8.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    18 October 2004
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Input, double X, the number whose sign is desired.
    //
    //    Output, double R8_SIGN, the sign of X.
    //
  {
    double value;

    if ( x < 0.0 )
    {
      value = -1.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }
  //****************************************************************************80

  void timestamp ( )

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    TIMESTAMP prints the current YMDHMS date as a time stamp.
    //
    //  Example:
    //
    //    31 May 2001 09:45:54 AM
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    24 September 2003
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    None
    //
  {
    const int TIME_SIZE(40);

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    time_t now;

    now = time ( NULL );
    tm = localtime ( &now );

    strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

    cout << time_buffer << "\n";

    return;
  }
  //****************************************************************************80

  // double local_min ( double a, double b, double t, func_base& f,
  //                    double &x )
  //
  //   //****************************************************************************80
  //   //
  //   //  Purpose:
  //   //
  //   //    LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
  //   //
  //   //  Discussion:
  //   //
  //   //    The method used is a combination of golden section search and
  //   //    successive parabolic interpolation.  Convergence is never much slower
  //   //    than that for a Fibonacci search.  If F has a continuous second
  //   //    derivative which is positive at the minimum (which is not at A or
  //   //    B), then convergence is superlinear, and usually of the order of
  //   //    about 1.324....
  //   //
  //   //    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
  //   //    F is never evaluated at two points closer than TOL.
  //   //
  //   //    If F is a unimodal function and the computed values of F are always
  //   //    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
  //   //    LOCAL_MIN approximates the abscissa of the global minimum of F on the
  //   //    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
  //   //
  //   //    If F is not unimodal, then LOCAL_MIN may approximate a local, but
  //   //    perhaps non-global, minimum to the same accuracy.
  //   //
  //   //  Licensing:
  //   //
  //   //    This code is distributed under the GNU LGPL license.
  //   //
  //   //  Modified:
  //   //
  //   //    17 July 2011
  //   //
  //   //  Author:
  //   //
  //   //    Original FORTRAN77 version by Richard Brent.
  //   //    C++ version by John Burkardt.
  //   //
  //   //  Reference:
  //   //
  //   //    Richard Brent,
  //   //    Algorithms for Minimization Without Derivatives,
  //   //    Dover, 2002,
  //   //    ISBN: 0-486-41998-3,
  //   //    LC: QA402.5.B74.
  //   //
  //   //  Parameters:
  //   //
  //   //    Input, double A, B, the endpoints of the interval.
  //   //
  //   //    Input, double T, a positive absolute error tolerance.
  //   //
  //   //    Input, func_base& F, a user-supplied c++ functor whose
  //   //    local minimum is being sought.  The input and output
  //   //    of F() are of type double.
  //   //
  //   //    Output, double &X, the estimated value of an abscissa
  //   //    for which F attains a local minimum value in [A,B].
  //   //
  //   //    Output, double LOCAL_MIN, the value F(X).
  //   //
  // {
  //   double c;
  //   double d;
  //   double e;
  //   double eps;
  //   double fu;
  //   double fv;
  //   double fw;
  //   double fx;
  //   double m;
  //   double p;
  //   double q;
  //   double r;
  //   double sa;
  //   double sb;
  //   double t2;
  //   double tol;
  //   double u;
  //   double v;
  //   double w;
  //   //
  //   //  C is the square of the inverse of the golden ratio.
  //   //
  //   c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );
  //
  //   eps = sqrt ( r8_epsilon ( ) );
  //
  //   sa = a;
  //   sb = b;
  //   x = sa + c * ( b - a );
  //   w = x;
  //   v = w;
  //   e = 0.0;
  //   fx = f ( x );
  //   fw = fx;
  //   fv = fw;
  //
  //   for ( ; ; )
  //   {
  //     m = 0.5 * ( sa + sb ) ;
  //     tol = eps * fabs ( x ) + t;
  //     t2 = 2.0 * tol;
  //     //
  //     //  Check the stopping criterion.
  //     //
  //     if ( fabs ( x - m ) <= t2 - 0.5 * ( sb - sa ) )
  //     {
  //       break;
  //     }
  //     //
  //     //  Fit a parabola.
  //     //
  //     r = 0.0;
  //     q = r;
  //     p = q;
  //
  //     if ( tol < fabs ( e ) )
  //     {
  //       r = ( x - w ) * ( fx - fv );
  //       q = ( x - v ) * ( fx - fw );
  //       p = ( x - v ) * q - ( x - w ) * r;
  //       q = 2.0 * ( q - r );
  //       if ( 0.0 < q )
  //       {
  //         p = - p;
  //       }
  //       q = fabs ( q );
  //       r = e;
  //       e = d;
  //     }
  //
  //     if ( fabs ( p ) < fabs ( 0.5 * q * r ) &&
  //          q * ( sa - x ) < p &&
  //          p < q * ( sb - x ) )
  //     {
  //       //
  //       //  Take the parabolic interpolation step.
  //       //
  //       d = p / q;
  //       u = x + d;
  //       //
  //       //  F must not be evaluated too close to A or B.
  //       //
  //       if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
  //       {
  //         if ( x < m )
  //         {
  //           d = tol;
  //         }
  //         else
  //         {
  //           d = - tol;
  //         }
  //       }
  //     }
  //     //
  //     //  A golden-section step.
  //     //
  //     else
  //     {
  //       if ( x < m )
  //       {
  //         e = sb - x;
  //       }
  //       else
  //       {
  //         e = sa - x;
  //       }
  //       d = c * e;
  //     }
  //     //
  //     //  F must not be evaluated too close to X.
  //     //
  //     if ( tol <= fabs ( d ) )
  //     {
  //       u = x + d;
  //     }
  //     else if ( 0.0 < d )
  //     {
  //       u = x + tol;
  //     }
  //     else
  //     {
  //       u = x - tol;
  //     }
  //
  //     fu = f ( u );
  //     //
  //     //  Update A, B, V, W, and X.
  //     //
  //     if ( fu <= fx )
  //     {
  //       if ( u < x )
  //       {
  //         sb = x;
  //       }
  //       else
  //       {
  //         sa = x;
  //       }
  //       v = w;
  //       fv = fw;
  //       w = x;
  //       fw = fx;
  //       x = u;
  //       fx = fu;
  //     }
  //     else
  //     {
  //       if ( u < x )
  //       {
  //         sa = u;
  //       }
  //       else
  //       {
  //         sb = u;
  //       }
  //
  //       if ( fu <= fw || w == x )
  //       {
  //         v = w;
  //         fv = fw;
  //         w = u;
  //         fw = fu;
  //       }
  //       else if ( fu <= fv || v == x || v == w )
  //       {
  //         v = u;
  //         fv = fu;
  //       }
  //     }
  //   }
  //   return fx;
  // }
  //

  //****************************************************************************80

  // double local_min ( double a, double b, double t, double f ( double x ),
  //   double &x ){
  //   func_wrapper foo(f);
  //   return local_min(a, b, t, foo, x);
  // }





} // end namespace brent



