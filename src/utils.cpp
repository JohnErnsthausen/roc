#include "utils.h"

#define x(i) x[i-1]
#define y(i) y[i-1]

int sgn(double val)
{
  return (int) ((double(0) < val) - (val < double(0)));
}

double sign(double a, double b)
{
  return (b>=0.0) ? abs(a) : -abs(a);
}

void swap(double *x, double *y)
{
  double tmp = *x;
  *x = *y;
  *y = tmp;
}

//  Generates an n-dimensional Householder reflector
//
//     H = I - tau*( 1 ) * ( 1 v' ),    H' * H = I,   tau scalar
//                 ( v ) 
//
//  such that
//
//     H * ( alpha ) = ( beta ),    alpha, beta scalars
//         (   x   )   (   0  )     x  (n-1)-dimensional vector
//
//  Because of H'* H = I it follows that
// 
//     alpha^2 + x^T x = beta^2    ==> beta = sqrt(alpha^2 + x^T x)
//
//     H'( beta )  = ( alpha )     ==> tau = (beta - alpha)/beta
//       (  0   )    (  x    )           v = xscal*x,  
//                                   xscal = 1/(beta - alpha)
//       
//     If x = 0, then tau = 0 and H = I, otherwise  1 <= tau <= 2.
//
//  This is an edited version of the LAPACK routine DLARFG
//
//  Variables in the calling sequence:
//  ----------------------------------
//  N      I   IN   Dimension of H
//  ALPHA  D   IN   The scalar alpha
//             OUT  The scalar beta
//  X      D   IN   The given vector x of dimension n - 1
//             OUT  The vector v
//  INCX   I   IN   The increment between elements of X, INCX .NE. 0
//  TAU    D   OUT  The scalar tau
//  SAFMIN D   IN   Safe minimum such that 1.0/SAFMIN does not
//                  overflow
int housg(int n, double *alpha, double *x, int incx, double *tau, double safmin)
{
  const double one=1.0, zer=0.0;
  int ix=0, knt=0, kx=0, nm1=0;
  double a1=0.0, a2=0.0, beta=0.0, rsafmn=0.0, tmp=0.0, xnorm=0.0, xscal=0.0;

  // H is defined as the identity whenever incx==0; implies don't iterate through x
  if( incx==0 )
  {
    *tau = zer;
    return -2;
  }

  // H is defined as the identity whenever n<1
  if( n<1 )
  {
    *tau = zer;
    return -1;
  }

  // H is defined as the identity whenever n=1
  if( n==1 )
  {
    *tau = zer;
    return 0;
  }

  // Dimension of x
  nm1 = n - 1;
  // Norm of x
  xnorm = ddist2(nm1, x, incx, x, incx, 1);

  // H is the identity whenever xnorm=0 with alpha=beta and tau=0
  if( xnorm==zer )
  {
    *tau = zer;
    return 0;
  }

  // General case
  kx = 1;
  if( incx<0 ) kx = 1 - (n-1)*incx;

  a1 = xnorm;
  a2 = abs(*alpha);
  if( a1<a2 )
  {
    a1 = a2;
    a2 = xnorm;
  }
  if( a2==zer )
  {
    beta = a1;
  }
  else
  {
    tmp = a1/a2;
    beta = a2*sqrt( one + tmp*tmp );
  }
  if( *alpha>zer ) beta = -beta;
  if( abs(beta)>=safmin )
  {
    *tau = (beta - *alpha) / beta;
    xscal = one / (*alpha - beta);
    *alpha = beta;
  }
  else // xnorm, beta may be inaccurate; scale x and recompute
  {
    rsafmn = one / safmin;
    knt = 0;
    do
    {
      knt = knt + 1;
      if( incx==1 )
      {
        for(int i=1;i<=nm1;i++)
        {
          x(i) = rsafmn*x(i);
        }
      }
      else
      {
        ix = kx;
        for(int i=1;i<=nm1;i++)
        {
          x(ix) = rsafmn*x(ix);
          ix = ix + incx;
        }
      }
      beta = beta*rsafmn;
      *alpha = (*alpha)*rsafmn;
    }
    while( abs( beta )<safmin );

    // new beta satisfies safmin <= beta <= 1.0
    xnorm = ddist2(nm1, x, incx, x, incx, 1);
    a1 = xnorm;
    a2 = abs(*alpha);
    if( a1<a2 )
    {
      a1 = a2;
      a2 = xnorm;
    }
    if( a2==zer )
    {
      beta = a1;
    }
    else
    {
      tmp = a1/a2;
      beta = a2*sqrt( one + tmp*tmp);
    }
    if( *alpha>zer ) beta = -beta;
    *tau = (beta - *alpha) / beta;
    xscal = one / (*alpha - beta);
    *alpha = beta;
    for(int j=1;j<=knt;j++)
    {
      *alpha = (*alpha)*safmin;
    }
  }
  if( incx==1 )
  {
    for(int i=1;i<=nm1;i++)
    {
      x(i) = xscal*x(i);
    }
  }
  else
  {
    ix = kx;
    for(int i=1;i<=nm1;i++)
    {
      x(ix) = xscal*x(ix);
      ix = ix + incx;
    }
  }
  return 0;
}

//  Computes either the Euclidean distance between two N-dimensional 
//  vectors X and Y or the Euclidean norm of one such vector X. 
//
//  Call the routine 
//  either
//     with KVEC = 2 and two vectors X and Y of dimension N stored
//     with storage-increments INCX and INCY, respectively.
//  or
//     with KVEC = 1 and one vector X of dimension N stored
//     with storage-increment INCX. In this case the second array 
//     is not referenced and can be a dummy array or simply the
//     array X again.
//
//  If N .LE. 0 then zero is returned, if N .GE. 1 then the
//  storage increments cannot be zero.
//
//  The algorithm follows the four-phase method of C. L. Lawson in the
//  LAPACK routine DNRM2.F.  As in DNRM2.F two built-in constants are 
//  used that are hopefully applicable to all machines.
//     CUTLO = maximum of  DSQRT(u/eps)  over all known machines.
//     CUTHI = minimum of  DSQRT(v)      over all known machines.
//  where
//     eps = smallest number such that 1.0d0 + eps .gt. 1.0d0
//     u   = smallest positive number  (underflow limit)
//     v   = largest  number           (overflow  limit)
//
//  Values for CUTLO and CUTHI listed in DNRM2.F are as follows:
//
//     CUTLO, s.p.  u/eps = 2**(-102) for honeywell. Close 
//                  seconds are univac and dec at 2**(-103)
//                  thus CUTLO = 2**(-51) = 4.44089e-16
//     CUTHI, s.p.  v = 2**127 for univac, honeywell, and dec.
//                  thus CUTHI = 2**(63.5) = 1.30438e19
//     CUTLO, d.p.  u/eps = 2**(-67) for honeywell and dec.
//                  thus CUTLO = 2**(-33.5) = 8.23181d-11
//     CUTHI, d.p.  same as s.p.  CUTHI = 1.30438d19
//
//     data cutlo, cuthi / 8.232d-11,  1.304d19 /
//     data cutlo, cuthi / 4.441e-16,  1.304e19 /                 
//
//  In line with the four phases of DNRM2.F, the algorithm uses 
//  four states identified by LEVEL = 0,1,2,3, respectively, 
//  which correspond to the following cases:
//
//     LEVEL = 0  only zero terms have been found so far
//     LEVEL = 1  all nonzero terms encountered so far do not
//                exceed CUTLO in modulus
//     LEVEL = 2  there are some terms that are larger than
//                CUTLO in modulus but none exceeds 
//                HITEST = CUTLO/DBLE(N) 
//     LEVEL = 3  there are terms that exceed HITEST in modulus.
//
//  All state transitions can only increase the LEVEL.
//
//  Variables in the calling sequence
//  ---------------------------------
//     N     I    IN   Dimension of the vectors X and Y
//     X     D    IN   The first vector of dimension N
//     INCX  I    IN   Storage increment of X
//     Y     D    IN   The second vector of dimension N
//     INYY  D    IN   Storage increment of Y
//     KVEC  I    IN   Number of vectors
//                     KVEC = 1  Only one vector, namely X, is given,
//                               the Euclidean norm of X is computed  
//                               and the Y array is not referenced
//                               This is the default
//                     KVEC = 2  Two vectors X and Y  are given,
//                               the Euclidean distance between 
//                               X and Y is computed
//
double ddist2( int n, double *x, int incx, double *y, int incy, int kvec)
{
  const double zer=0.0, one=1.0, cutlo=8.232e-11, cuthi=1.304e19;

  int ix=1, iy=1, level=0, job=1;
  double ddist2 = zer;
  double diff=0.0, dx=0.0, dy=0.0;
  double hitest=0.0;
  double sum=0.0, tmp=0.0, trm=0.0, xmax=0.0;

  // The inner product is computed whenever kvec=2.
  // Otherwise proceeds as if kvec=1.
  job= ( kvec==2 ) ? 2 : 1;

  // If n<=0, then define the norm of x or the inner product of x with y to be zero
  if( n <= 0 ) return zer;

  // If the storage-increment incx equals zero, then define the norm of x to be zero
  if( job==1 && incx==0 ) return zer;

  // If both incx and incy equal zero, then define the inner product of x with y to be zero
  if( job==2 && incx==0 && incy==0 ) return zer;
  
  // Initializations

  hitest = cuthi/n;
  if( hitest < cutlo ) hitest = cutlo;

  ix = 1;
  if( incx < 0 ) ix = (-n+1)*incx + 1;
  iy = 1;
  if( incy < 0 ) iy = (-n+1)*incy + 1;

  sum = zer;
  xmax = zer;
  level = 0;

  for(int j=1;j<n;j++)
  {
    dx = x(ix);
    ix = ix + incx;
    if( job==1 )  // One vector case
    {
      diff = dx;
    }
    else  // Two vector case
    {
      dy = y(iy);
      iy = iy + incy;

      // Compute diff = dx - dy to avoid subtractive cancellation

      if( sgn(dx) != sgn(dy) )
      {
        // If dx and dy are not of the same sgn, then safely subtract
        diff = dx - dy;
      }
      else
      {
        // If dx and dy are of the same sign, then subtract by
        // 1.) make sure dx > dy
        // 2.a) We have sgn(dx)=sgn(dy).
        //      If dx=0, then dy=0. Thus diff=0.
        // 2.b) else Subtract dx - dy = dx*(one-dy/dx)
        if( abs(dx) < abs(dy) )
        {
           tmp = dy;
           dy  = dx;
           dx  = tmp;
        }
        if( dx == zer )
        {
           diff = zer;
        }
        else
        {
           diff = dx*(one - dy/dx);
        }
      }
    }
    
    // Sum the squares of the nonzero terms

    if( diff!=zer )
    {
      trm = abs(diff);
      if( trm<=cutlo ) // Very small terms
      {
        if( level==0 ) // The first nonzero term encountered
        {
          level = 1;
          xmax = trm;
          tmp = trm/xmax;
          sum = sum + tmp*tmp;
        }
        else
        { 
          if( level==1 )
          {
            if( trm>xmax )
            {
              tmp = xmax/trm;
              sum = one + sum*tmp*tmp;
              xmax = trm;
            }
            else
            {
              tmp = trm/xmax;
              sum = sum + tmp*tmp;
            }
          }
          else
          {           
            sum = sum + trm*trm;
          }
        }
      }
      else
      {
        // Mid-sized terms -- transition to level 2
        if( level==0 )
        {
           level = 2;
        }
        else
        {
          if( level==1)
          {
            level = 2;
            sum = (sum * xmax) * xmax;
          }
        }
  
        if( trm<=hitest )
        {
          sum = sum + trm*trm;
        }
        else
        {
          // large terms
          if( level==2 )
          {
            // transition to level = 3
            level = 3;
            sum = (sum / trm) / trm;
            xmax = trm;
          }

          if( trm>xmax )
          {
            tmp = xmax/trm;
            sum = one + sum*tmp*tmp;
            xmax = trm;
          }
          else
          {
            tmp = trm/xmax;
            sum = sum + tmp*tmp;
          }
        }
      } // Mid
    }
  }

  if( (xmax==zer) || (level==2) )
  {
    ddist2 = sqrt(sum);
  }
  else
  {
    ddist2 = xmax *sqrt(sum);
  }
  return ddist2;
}
