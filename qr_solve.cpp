# include <cmath>
# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "qr_solve.hpp"
# include "r8lib.hpp"

//****************************************************************************80

void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DAXPY computes constant times a vector plus a vector.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in DX and DY.
//
//    Input, double DA, the multiplier of DX.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries of DX.
//
//    Input/output, double DY[*], the second vector.
//    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
//
//    Input, int INCY, the increment between successive entries of DY.
//
{
  int i;
  int ix;
  int iy;
  int m;

  if ( n <= 0 )
  {
    return;
  }

  if ( da == 0.0 )
  {
    return;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dy[iy] = dy[iy] + da * dx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 4;

    for ( i = 0; i < m; i++ )
    {
      dy[i] = dy[i] + da * dx[i];
    }

    for ( i = m; i < n; i = i + 4 )
    {
      dy[i  ] = dy[i  ] + da * dx[i  ];
      dy[i+1] = dy[i+1] + da * dx[i+1];
      dy[i+2] = dy[i+2] + da * dx[i+2];
      dy[i+3] = dy[i+3] + da * dx[i+3];
    }

  }

  return;
}
//****************************************************************************80

double ddot ( int n, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DDOT forms the dot product of two vectors.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries in DX.
//
//    Input, double DY[*], the second vector.
//
//    Input, int INCY, the increment between successive entries in DY.
//
//    Output, double DDOT, the sum of the product of the corresponding
//    entries of DX and DY.
//
{
  double dtemp;
  int i;
  int ix;
  int iy;
  int m;

  dtemp = 0.0;

  if ( n <= 0 )
  {
    return dtemp;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dtemp = dtemp + dx[ix] * dy[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      dtemp = dtemp + dx[i] * dy[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      dtemp = dtemp + dx[i  ] * dy[i  ]
                    + dx[i+1] * dy[i+1]
                    + dx[i+2] * dy[i+2]
                    + dx[i+3] * dy[i+3]
                    + dx[i+4] * dy[i+4];
    }

  }

  return dtemp;
}
//****************************************************************************80

double dnrm2 ( int n, double x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    DNRM2 returns the euclidean norm of a vector.
//
//  Discussion:
//
//     DNRM2 ( X ) = sqrt ( X' * X )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[*], the vector whose norm is to be computed.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Output, double DNRM2, the Euclidean norm of X.
//
{
  double absxi;
  int i;
  int ix;
  double norm;
  double scale;
  double ssq;

  if ( n < 1 || incx < 1 )
  {
    norm = 0.0;
  }
  else if ( n == 1 )
  {
    norm = fabs ( x[0] );
  }
  else
  {
    scale = 0.0;
    ssq = 1.0;
    ix = 0;

    for ( i = 0; i < n; i++ )
    {
      if ( x[ix] != 0.0 )
      {
        absxi = fabs ( x[ix] );
        if ( scale < absxi )
        {
          ssq = 1.0 + ssq * ( scale / absxi ) * ( scale / absxi );
          scale = absxi;
        }
        else
        {
          ssq = ssq + ( absxi / scale ) * ( absxi / scale );
        }
      }
      ix = ix + incx;
    }

    norm  = scale * sqrt ( ssq );
  }

  return norm;
}
//****************************************************************************80

void dqrank ( double a[], int lda, int m, int n, double tol, int &kr, 
  int jpvt[], double qraux[] )

//****************************************************************************80
//
//  Purpose:
//
//    DQRANK computes the QR factorization of a rectangular matrix.
//
//  Discussion:
//
//    This routine is used in conjunction with sqrlss to solve
//    overdetermined, underdetermined and singular linear systems
//    in a least squares sense.
//
//    DQRANK uses the LINPACK subroutine DQRDC to compute the QR
//    factorization, with column pivoting, of an M by N matrix A.
//    The numerical rank is determined using the tolerance TOL.
//
//    Note that on output, ABS ( A(1,1) ) / ABS ( A(KR,KR) ) is an estimate
//    of the condition number of the matrix of independent columns,
//    and of R.  This estimate will be <= 1/TOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the matrix whose
//    decomposition is to be computed.  On output, the information from DQRDC.
//    The triangular matrix R of the QR factorization is contained in the
//    upper triangle and information needed to recover the orthogonal
//    matrix Q is stored below the diagonal in A and in the vector QRAUX.
//
//    Input, int LDA, the leading dimension of A, which must
//    be at least M.
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input, double TOL, a relative tolerance used to determine the
//    numerical rank.  The problem should be scaled so that all the elements
//    of A have roughly the same absolute accuracy, EPS.  Then a reasonable
//    value for TOL is roughly EPS divided by the magnitude of the largest
//    element.
//
//    Output, int &KR, the numerical rank.
//
//    Output, int JPVT[N], the pivot information from DQRDC.
//    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
//    independent to within the tolerance TOL and the remaining columns
//    are linearly dependent.
//
//    Output, double QRAUX[N], will contain extra information defining
//    the QR factorization.
//
{
  int i;
  int j;
  int job;
  int k;
  double *work;

  for ( i = 0; i < n; i++ )
  {
    jpvt[i] = 0;
  }

  work = new double[n];
  job = 1;

  dqrdc ( a, lda, m, n, qraux, jpvt, work, job );

  kr = 0;
  k = i4_min ( m, n );

  for ( j = 0; j < k; j++ )
  {
    if ( fabs ( a[j+j*lda] ) <= tol * fabs ( a[0+0*lda] ) )
    {
      return;
    }
    kr = j + 1;
  }

  delete [] work;

  return;
}
//****************************************************************************80

void dqrdc ( double a[], int lda, int n, int p, double qraux[], int jpvt[], 
  double work[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DQRDC computes the QR factorization of a real rectangular matrix.
//
//  Discussion:
//
//    DQRDC uses Householder transformations.
//
//    Column pivoting based on the 2-norms of the reduced columns may be
//    performed at the user's option.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2005
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A(LDA,P).  On input, the N by P matrix
//    whose decomposition is to be computed.  On output, A contains in
//    its upper triangle the upper triangular matrix R of the QR
//    factorization.  Below its diagonal A contains information from
//    which the orthogonal part of the decomposition can be recovered.
//    Note that if pivoting has been requested, the decomposition is not that
//    of the original matrix A but that of A with its columns permuted
//    as described by JPVT.
//
//    Input, int LDA, the leading dimension of the array A.  LDA must
//    be at least N.
//
//    Input, int N, the number of rows of the matrix A.
//
//    Input, int P, the number of columns of the matrix A.
//
//    Output, double QRAUX[P], contains further information required
//    to recover the orthogonal part of the decomposition.
//
//    Input/output, integer JPVT[P].  On input, JPVT contains integers that
//    control the selection of the pivot columns.  The K-th column A(*,K) of A
//    is placed in one of three classes according to the value of JPVT(K).
//      > 0, then A(K) is an initial column.
//      = 0, then A(K) is a free column.
//      < 0, then A(K) is a final column.
//    Before the decomposition is computed, initial columns are moved to
//    the beginning of the array A and final columns to the end.  Both
//    initial and final columns are frozen in place during the computation
//    and only free columns are moved.  At the K-th stage of the
//    reduction, if A(*,K) is occupied by a free column it is interchanged
//    with the free column of largest reduced norm.  JPVT is not referenced
//    if JOB == 0.  On output, JPVT(K) contains the index of the column of the
//    original matrix that has been interchanged into the K-th column, if
//    pivoting was requested.
//
//    Workspace, double WORK[P].  WORK is not referenced if JOB == 0.
//
//    Input, int JOB, initiates column pivoting.
//    0, no pivoting is done.
//    nonzero, pivoting is done.
//
{
  int j;
  int jp;
  int l;
  int lup;
  int maxj;
  double maxnrm;
  double nrmxl;
  int pl;
  int pu;
  bool swapj;
  double t;
  double tt;

  pl = 1;
  pu = 0;
//
//  If pivoting is requested, rearrange the columns.
//
  if ( job != 0 )
  {
    for ( j = 1; j <= p; j++ )
    {
      swapj = ( 0 < jpvt[j-1] );

      if ( jpvt[j-1] < 0 )
      {
        jpvt[j-1] = -j;
      }
      else
      {
        jpvt[j-1] = j;
      }

      if ( swapj )
      {
        if ( j != pl )
        {
          dswap ( n, a+0+(pl-1)*lda, 1, a+0+(j-1), 1 );
        }
        jpvt[j-1] = jpvt[pl-1];
        jpvt[pl-1] = j;
        pl = pl + 1;
      }
    }
    pu = p;

    for ( j = p; 1 <= j; j-- )
    {
      if ( jpvt[j-1] < 0 )
      {
        jpvt[j-1] = -jpvt[j-1];

        if ( j != pu )
        {
          dswap ( n, a+0+(pu-1)*lda, 1, a+0+(j-1)*lda, 1 );
          jp = jpvt[pu-1];
          jpvt[pu-1] = jpvt[j-1];
          jpvt[j-1] = jp;
        }
        pu = pu - 1;
      }
    }
  }
//
//  Compute the norms of the free columns.
//
  for ( j = pl; j <= pu; j++ )
  {
    qraux[j-1] = dnrm2 ( n, a+0+(j-1)*lda, 1 );
  }

  for ( j = pl; j <= pu; j++ )
  {
    work[j-1] = qraux[j-1];
  }
//
//  Perform the Householder reduction of A.
//
  lup = i4_min ( n, p );

  for ( l = 1; l <= lup; l++ )
  {
//
//  Bring the column of largest norm into the pivot position.
//
    if ( pl <= l && l < pu )
    {
      maxnrm = 0.0;
      maxj = l;
      for ( j = l; j <= pu; j++ )
      {
        if ( maxnrm < qraux[j-1] )
        {
          maxnrm = qraux[j-1];
          maxj = j;
        }
      }

      if ( maxj != l )
      {
        dswap ( n, a+0+(l-1)*lda, 1, a+0+(maxj-1)*lda, 1 );
        qraux[maxj-1] = qraux[l-1];
        work[maxj-1] = work[l-1];
        jp = jpvt[maxj-1];
        jpvt[maxj-1] = jpvt[l-1];
        jpvt[l-1] = jp;
      }
    }
//
//  Compute the Householder transformation for column L.
//
    qraux[l-1] = 0.0;

    if ( l != n )
    {
      nrmxl = dnrm2 ( n-l+1, a+l-1+(l-1)*lda, 1 );

      if ( nrmxl != 0.0 )
      {
        if ( a[l-1+(l-1)*lda] != 0.0 )
        {
          nrmxl = nrmxl * r8_sign ( a[l-1+(l-1)*lda] );
        }

        dscal ( n-l+1, 1.0 / nrmxl, a+l-1+(l-1)*lda, 1 );
        a[l-1+(l-1)*lda] = 1.0 + a[l-1+(l-1)*lda];
//
//  Apply the transformation to the remaining columns, updating the norms.
//
        for ( j = l + 1; j <= p; j++ )
        {
          t = -ddot ( n-l+1, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 ) 
            / a[l-1+(l-1)*lda];
          daxpy ( n-l+1, t, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 );

          if ( pl <= j && j <= pu )
          {
            if ( qraux[j-1] != 0.0 )
            {
              tt = 1.0 - pow ( fabs ( a[l-1+(j-1)*lda] ) / qraux[j-1], 2 );
              tt = r8_max ( tt, 0.0 );
              t = tt;
              tt = 1.0 + 0.05 * tt * pow ( qraux[j-1] / work[j-1], 2 );

              if ( tt != 1.0 )
              {
                qraux[j-1] = qraux[j-1] * sqrt ( t );
              }
              else
              {
                qraux[j-1] = dnrm2 ( n-l, a+l+(j-1)*lda, 1 );
                work[j-1] = qraux[j-1];
              }
            }
          }
        }
//
//  Save the transformation.
//
        qraux[l-1] = a[l-1+(l-1)*lda];
        a[l-1+(l-1)*lda] = -nrmxl;
      }
    }
  }
  return;
}
//****************************************************************************80

int dqrls ( double a[], int lda, int m, int n, double tol, int &kr, double b[], 
  double x[], double rsd[], int jpvt[], double qraux[], int itask )

//****************************************************************************80
//
//  Purpose:
//
//    DQRLS factors and solves a linear system in the least squares sense.
//
//  Discussion:
//
//    The linear system may be overdetermined, underdetermined or singular.
//    The solution is obtained using a QR factorization of the
//    coefficient matrix.
//
//    DQRLS can be efficiently used to solve several least squares
//    problems with the same matrix A.  The first system is solved
//    with ITASK = 1.  The subsequent systems are solved with
//    ITASK = 2, to avoid the recomputation of the matrix factors.
//    The parameters KR, JPVT, and QRAUX must not be modified
//    between calls to DQRLS.
//
//    DQRLS is used to solve in a least squares sense
//    overdetermined, underdetermined and singular linear systems.
//    The system is A*X approximates B where A is M by N.
//    B is a given M-vector, and X is the N-vector to be computed.
//    A solution X is found which minimimzes the sum of squares (2-norm)
//    of the residual,  A*X - B.
//
//    The numerical rank of A is determined using the tolerance TOL.
//
//    DQRLS uses the LINPACK subroutine DQRDC to compute the QR
//    factorization, with column pivoting, of an M by N matrix A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters:
//
//    Input/output, double A[LDA*N], an M by N matrix.
//    On input, the matrix whose decomposition is to be computed.
//    In a least squares data fitting problem, A(I,J) is the
//    value of the J-th basis (model) function at the I-th data point.
//    On output, A contains the output from DQRDC.  The triangular matrix R
//    of the QR factorization is contained in the upper triangle and
//    information needed to recover the orthogonal matrix Q is stored
//    below the diagonal in A and in the vector QRAUX.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input, double TOL, a relative tolerance used to determine the
//    numerical rank.  The problem should be scaled so that all the elements
//    of A have roughly the same absolute accuracy EPS.  Then a reasonable
//    value for TOL is roughly EPS divided by the magnitude of the largest
//    element.
//
//    Output, int &KR, the numerical rank.
//
//    Input, double B[M], the right hand side of the linear system.
//
//    Output, double X[N], a least squares solution to the linear
//    system.
//
//    Output, double RSD[M], the residual, B - A*X.  RSD may
//    overwrite B.
//
//    Workspace, int JPVT[N], required if ITASK = 1.
//    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
//    independent to within the tolerance TOL and the remaining columns
//    are linearly dependent.  ABS ( A(1,1) ) / ABS ( A(KR,KR) ) is an estimate
//    of the condition number of the matrix of independent columns,
//    and of R.  This estimate will be <= 1/TOL.
//
//    Workspace, double QRAUX[N], required if ITASK = 1.
//
//    Input, int ITASK.
//    1, DQRLS factors the matrix A and solves the least squares problem.
//    2, DQRLS assumes that the matrix A was factored with an earlier
//       call to DQRLS, and only solves the least squares problem.
//
//    Output, int DQRLS, error code.
//    0:  no error
//    -1: LDA < M   (fatal error)
//    -2: N < 1     (fatal error)
//    -3: ITASK < 1 (fatal error)
//
{
  int ind;

  if ( lda < m )
  {
    cerr << "\n";
    cerr << "DQRLS - Fatal error!\n";
    cerr << "  LDA < M.\n";
    ind = -1;
    return ind;
  }

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "DQRLS - Fatal error!\n";
    cerr << "  N <= 0.\n";
    ind = -2;
    return ind;
  }

  if ( itask < 1 )
  {
    cerr << "\n";
    cerr << "DQRLS - Fatal error!\n";
    cerr << "  ITASK < 1.\n";
    ind = -3;
    return ind;
  }

  ind = 0;
//
//  Factor the matrix.
//
  if ( itask == 1 )
  {
    dqrank ( a, lda, m, n, tol, kr, jpvt, qraux );
  }
//
//  Solve the least-squares problem.
//
  dqrlss ( a, lda, m, n, kr, b, x, rsd, jpvt, qraux );

  return ind;
}
//****************************************************************************80

void dqrlss ( double a[], int lda, int m, int n, int kr, double b[], double x[], 
  double rsd[], int jpvt[], double qraux[] )

//****************************************************************************80
//
//  Purpose:
//
//    DQRLSS solves a linear system in a least squares sense.
//
//  Discussion:
//
//    DQRLSS must be preceeded by a call to DQRANK.
//
//    The system is to be solved is
//      A * X = B
//    where
//      A is an M by N matrix with rank KR, as determined by DQRANK,
//      B is a given M-vector,
//      X is the N-vector to be computed.
//
//    A solution X, with at most KR nonzero components, is found which
//    minimizes the 2-norm of the residual (A*X-B).
//
//    Once the matrix A has been formed, DQRANK should be
//    called once to decompose it.  Then, for each right hand
//    side B, DQRLSS should be called once to obtain the
//    solution and residual.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2012
//
//  Parameters:
//
//    Input, double A[LDA*N], the QR factorization information
//    from DQRANK.  The triangular matrix R of the QR factorization is
//    contained in the upper triangle and information needed to recover
//    the orthogonal matrix Q is stored below the diagonal in A and in
//    the vector QRAUX.
//
//    Input, int LDA, the leading dimension of A, which must
//    be at least M.
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input, int KR, the rank of the matrix, as estimated by DQRANK.
//
//    Input, double B[M], the right hand side of the linear system.
//
//    Output, double X[N], a least squares solution to the
//    linear system.
//
//    Output, double RSD[M], the residual, B - A*X.  RSD may
//    overwite B.
//
//    Input, int JPVT[N], the pivot information from DQRANK.
//    Columns JPVT[0], ..., JPVT[KR-1] of the original matrix are linearly
//    independent to within the tolerance TOL and the remaining columns
//    are linearly dependent.
//
//    Input, double QRAUX[N], auxiliary information from DQRANK
//    defining the QR factorization.
//
{
  int i;
  int j;
  int job;
  int k;
  double t;

  if ( kr != 0 )
  {
    job = 110;
    dqrsl ( a, lda, m, kr, qraux, b, rsd, rsd, x, rsd, rsd, job );
  }

  for ( i = 0; i < n; i++ )
  {
    jpvt[i] = - jpvt[i];
  }

  for ( i = kr; i < n; i++ )
  {
    x[i] = 0.0;
  }

  for ( j = 1; j <= n; j++ )
  {
    if ( jpvt[j-1] <= 0 )
    {
      k = - jpvt[j-1];
      jpvt[j-1] = k;

      while ( k != j )
      {
        t = x[j-1];
        x[j-1] = x[k-1];
        x[k-1] = t;
        jpvt[k-1] = -jpvt[k-1];
        k = jpvt[k-1];
      }
    }
  }
  return;
}
//****************************************************************************80

int dqrsl ( double a[], int lda, int n, int k, double qraux[], double y[], 
  double qy[], double qty[], double b[], double rsd[], double ab[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DQRSL computes transformations, projections, and least squares solutions.
//
//  Discussion:
//
//    DQRSL requires the output of DQRDC.
//
//    For K <= min(N,P), let AK be the matrix
//
//      AK = ( A(JPVT[0]), A(JPVT(2)), ..., A(JPVT(K)) )
//
//    formed from columns JPVT[0], ..., JPVT(K) of the original
//    N by P matrix A that was input to DQRDC.  If no pivoting was
//    done, AK consists of the first K columns of A in their
//    original order.  DQRDC produces a factored orthogonal matrix Q
//    and an upper triangular matrix R such that
//
//      AK = Q * (R)
//               (0)
//
//    This information is contained in coded form in the arrays
//    A and QRAUX.
//
//    The parameters QY, QTY, B, RSD, and AB are not referenced
//    if their computation is not requested and in this case
//    can be replaced by dummy variables in the calling program.
//    To save storage, the user may in some cases use the same
//    array for different parameters in the calling sequence.  A
//    frequently occuring example is when one wishes to compute
//    any of B, RSD, or AB and does not need Y or QTY.  In this
//    case one may identify Y, QTY, and one of B, RSD, or AB, while
//    providing separate arrays for anything else that is to be
//    computed.
//
//    Thus the calling sequence
//
//      dqrsl ( a, lda, n, k, qraux, y, dum, y, b, y, dum, 110, info )
//
//    will result in the computation of B and RSD, with RSD
//    overwriting Y.  More generally, each item in the following
//    list contains groups of permissible identifications for
//    a single calling sequence.
//
//      1. (Y,QTY,B) (RSD) (AB) (QY)
//
//      2. (Y,QTY,RSD) (B) (AB) (QY)
//
//      3. (Y,QTY,AB) (B) (RSD) (QY)
//
//      4. (Y,QY) (QTY,B) (RSD) (AB)
//
//      5. (Y,QY) (QTY,RSD) (B) (AB)
//
//      6. (Y,QY) (QTY,AB) (B) (RSD)
//
//    In any group the value returned in the array allocated to
//    the group corresponds to the last member of the group.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2005
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double A[LDA*P], contains the output of DQRDC.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the number of rows of the matrix AK.  It must
//    have the same value as N in DQRDC.
//
//    Input, int K, the number of columns of the matrix AK.  K
//    must not be greater than min(N,P), where P is the same as in the
//    calling sequence to DQRDC.
//
//    Input, double QRAUX[P], the auxiliary output from DQRDC.
//
//    Input, double Y[N], a vector to be manipulated by DQRSL.
//
//    Output, double QY[N], contains Q * Y, if requested.
//
//    Output, double QTY[N], contains Q' * Y, if requested.
//
//    Output, double B[K], the solution of the least squares problem
//      minimize norm2 ( Y - AK * B),
//    if its computation has been requested.  Note that if pivoting was
//    requested in DQRDC, the J-th component of B will be associated with
//    column JPVT(J) of the original matrix A that was input into DQRDC.
//
//    Output, double RSD[N], the least squares residual Y - AK * B,
//    if its computation has been requested.  RSD is also the orthogonal
//    projection of Y onto the orthogonal complement of the column space
//    of AK.
//
//    Output, double AB[N], the least squares approximation Ak * B,
//    if its computation has been requested.  AB is also the orthogonal
//    projection of Y onto the column space of A.
//
//    Input, integer JOB, specifies what is to be computed.  JOB has
//    the decimal expansion ABCDE, with the following meaning:
//
//      if A != 0, compute QY.
//      if B != 0, compute QTY.
//      if C != 0, compute QTY and B.
//      if D != 0, compute QTY and RSD.
//      if E != 0, compute QTY and AB.
//
//    Note that a request to compute B, RSD, or AB automatically triggers
//    the computation of QTY, for which an array must be provided in the
//    calling sequence.
//
//    Output, int DQRSL, is zero unless the computation of B has
//    been requested and R is exactly singular.  In this case, INFO is the
//    index of the first zero diagonal element of R, and B is left unaltered.
//
{
  bool cab;
  bool cb;
  bool cqty;
  bool cqy;
  bool cr;
  int i;
  int info;
  int j;
  int jj;
  int ju;
  double t;
  double temp;
//
//  set info flag.
//
  info = 0;
//
//  Determine what is to be computed.
//
  cqy =  (   job / 10000          != 0 );
  cqty = ( ( job %  10000 )       != 0 );
  cb =   ( ( job %   1000 ) / 100 != 0 );
  cr =   ( ( job %    100 ) /  10 != 0 );
  cab =  ( ( job %     10 )       != 0 );

  ju = i4_min ( k, n-1 );
//
//  Special action when N = 1.
//
  if ( ju == 0 )
  {
    if ( cqy )
    {
      qy[0] = y[0];
    }

    if ( cqty )
    {
      qty[0] = y[0];
    }

    if ( cab )
    {
      ab[0] = y[0];
    }

    if ( cb )
    {
      if ( a[0+0*lda] == 0.0 )
      {
        info = 1;
      }
      else
      {
        b[0] = y[0] / a[0+0*lda];
      }
    }

    if ( cr )
    {
      rsd[0] = 0.0;
    }
    return info;
  }
//
//  Set up to compute QY or QTY.
//
  if ( cqy )
  {
    for ( i = 1; i <= n; i++ )
    {
      qy[i-1] = y[i-1];
    }
  }

  if ( cqty )
  {
    for ( i = 1; i <= n; i++ )
    {
      qty[i-1] = y[i-1];
    }
  }
//
//  Compute QY.
//
  if ( cqy )
  {
    for ( jj = 1; jj <= ju; jj++ )
    {
      j = ju - jj + 1;

      if ( qraux[j-1] != 0.0 )
      {
        temp = a[j-1+(j-1)*lda];
        a[j-1+(j-1)*lda] = qraux[j-1];
        t = -ddot ( n-j+1, a+j-1+(j-1)*lda, 1, qy+j-1, 1 ) / a[j-1+(j-1)*lda];
        daxpy ( n-j+1, t, a+j-1+(j-1)*lda, 1, qy+j-1, 1 );
        a[j-1+(j-1)*lda] = temp;
      }
    }
  }
//
//  Compute Q'*Y.
//
  if ( cqty )
  {
    for ( j = 1; j <= ju; j++ )
    {
      if ( qraux[j-1] != 0.0 )
      {
        temp = a[j-1+(j-1)*lda];
        a[j-1+(j-1)*lda] = qraux[j-1];
        t = -ddot ( n-j+1, a+j-1+(j-1)*lda, 1, qty+j-1, 1 ) / a[j-1+(j-1)*lda];
        daxpy ( n-j+1, t, a+j-1+(j-1)*lda, 1, qty+j-1, 1 );
        a[j-1+(j-1)*lda] = temp;
      }
    }
  }
//
//  Set up to compute B, RSD, or AB.
//
  if ( cb )
  {
    for ( i = 1; i <= k; i++ )
    {
      b[i-1] = qty[i-1];
    }
  }

  if ( cab )
  {
    for ( i = 1; i <= k; i++ )
    {
      ab[i-1] = qty[i-1];
    }
  }

  if ( cr && k < n )
  {
    for ( i = k+1; i <= n; i++ )
    {
      rsd[i-1] = qty[i-1];
    }
  }

  if ( cab && k+1 <= n )
  {
    for ( i = k+1; i <= n; i++ )
    {
      ab[i-1] = 0.0;
    }
  }

  if ( cr )
  {
    for ( i = 1; i <= k; i++ )
    {
      rsd[i-1] = 0.0;
    }
  }
//
//  Compute B.
//
  if ( cb )
  {
    for ( jj = 1; jj <= k; jj++ )
    {
      j = k - jj + 1;

      if ( a[j-1+(j-1)*lda] == 0.0 )
      {
        info = j;
        break;
      }

      b[j-1] = b[j-1] / a[j-1+(j-1)*lda];

      if ( j != 1 )
      {
        t = -b[j-1];
        daxpy ( j-1, t, a+0+(j-1)*lda, 1, b, 1 );
      }
    }
  }
//
//  Compute RSD or AB as required.
//
  if ( cr || cab )
  {
    for ( jj = 1; jj <= ju; jj++ )
    {
      j = ju - jj + 1;

      if ( qraux[j-1] != 0.0 )
      {
        temp = a[j-1+(j-1)*lda];
        a[j-1+(j-1)*lda] = qraux[j-1];

        if ( cr )
        {
          t = -ddot ( n-j+1, a+j-1+(j-1)*lda, 1, rsd+j-1, 1 ) 
            / a[j-1+(j-1)*lda];
          daxpy ( n-j+1, t, a+j-1+(j-1)*lda, 1, rsd+j-1, 1 );
        }

        if ( cab )
        {
          t = -ddot ( n-j+1, a+j-1+(j-1)*lda, 1, ab+j-1, 1 ) 
            / a[j-1+(j-1)*lda];
          daxpy ( n-j+1, t, a+j-1+(j-1)*lda, 1, ab+j-1, 1 );
        }
        a[j-1+(j-1)*lda] = temp;
      }
    }
  }

  return info;
}
//****************************************************************************80

void drot ( int n, double x[], int incx, double y[], int incy, double c,
  double s )

//****************************************************************************80
//
//  Purpose:
//
//    DROT applies a plane rotation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double X[*], one of the vectors to be rotated.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], one of the vectors to be rotated.
//
//    Input, int INCY, the increment between successive elements of Y.
//
//    Input, double C, S, parameters (presumably the cosine and
//    sine of some angle) that define a plane rotation.
//
{
  int i;
  int ix;
  int iy;
  double stemp;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      stemp = c * x[i] + s * y[i];
      y[i]  = c * y[i] - s * x[i];
      x[i]  = stemp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      stemp = c * x[ix] + s * y[iy];
      y[iy] = c * y[iy] - s * x[ix];
      x[ix] = stemp;
      ix = ix + incx;
      iy = iy + incy;
    }

  }

  return;
}
//****************************************************************************80

void drotg ( double *sa, double *sb, double *c, double *s )

//****************************************************************************80
//
//  Purpose:
//
//    DROTG constructs a Givens plane rotation.
//
//  Discussion:
//
//    Given values A and B, this routine computes
//
//    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
//          = sign ( B ) if abs ( A ) <= abs ( B );
//
//    R     = SIGMA * ( A * A + B * B );
//
//    C = A / R if R is not 0
//      = 1     if R is 0;
//
//    S = B / R if R is not 0,
//        0     if R is 0.
//
//    The computed numbers then satisfy the equation
//
//    (  C  S ) ( A ) = ( R )
//    ( -S  C ) ( B ) = ( 0 )
//
//    The routine also computes
//
//    Z = S     if abs ( A ) > abs ( B ),
//      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
//      = 1     if C is 0.
//
//    The single value Z encodes C and S, and hence the rotation:
//
//    If Z = 1, set C = 0 and S = 1;
//    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
//    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input/output, double *SA, *SB,  On input, SA and SB are the values
//    A and B.  On output, SA is overwritten with R, and SB is
//    overwritten with Z.
//
//    Output, double *C, *S, the cosine and sine of the Givens rotation.
//
{
  double r;
  double roe;
  double scale;
  double z;

  if ( fabs ( *sb ) < fabs ( *sa ) )
  {
    roe = *sa;
  }
  else
  {
    roe = *sb;
  }

  scale = fabs ( *sa ) + fabs ( *sb );

  if ( scale == 0.0 )
  {
    *c = 1.0;
    *s = 0.0;
    r = 0.0;
  }
  else
  {
    r = scale * sqrt ( ( *sa / scale ) * ( *sa / scale )
                     + ( *sb / scale ) * ( *sb / scale ) );
    r = r8_sign ( roe ) * r;
    *c = *sa / r;
    *s = *sb / r;
  }

  if ( 0.0 < fabs ( *c ) && fabs ( *c ) <= *s )
  {
    z = 1.0 / *c;
  }
  else
  {
    z = *s;
  }

  *sa = r;
  *sb = z;

  return;
}
//****************************************************************************80

void dscal ( int n, double sa, double x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    DSCAL scales a vector by a constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double SA, the multiplier.
//
//    Input/output, double X[*], the vector to be scaled.
//
//    Input, int INCX, the increment between successive entries of X.
//
{
  int i;
  int ix;
  int m;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 )
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      x[i] = sa * x[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      x[i]   = sa * x[i];
      x[i+1] = sa * x[i+1];
      x[i+2] = sa * x[i+2];
      x[i+3] = sa * x[i+3];
      x[i+4] = sa * x[i+4];
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    for ( i = 0; i < n; i++ )
    {
      x[ix] = sa * x[ix];
      ix = ix + incx;
    }

  }

  return;
}
//****************************************************************************80

int dsvdc ( double a[], int lda, int m, int n, double s[], double e[], 
  double u[], int ldu, double v[], int ldv, double work[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DSVDC computes the singular value decomposition of a real rectangular matrix.
//
//  Discussion:
//
//    This routine reduces an M by N matrix A to diagonal form by orthogonal
//    transformations U and V.  The diagonal elements S(I) are the singular
//    values of A.  The columns of U are the corresponding left singular
//    vectors, and the columns of V the right singular vectors.
//
//    The form of the singular value decomposition is then
//
//      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2007
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the M by N matrix whose
//    singular value decomposition is to be computed.  On output, the matrix
//    has been destroyed.  Depending on the user's requests, the matrix may 
//    contain other useful information.
//
//    Input, int LDA, the leading dimension of the array A.
//    LDA must be at least M.
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix A.
//
//    Output, double S[MM], where MM = min(M+1,N).  The first
//    min(M,N) entries of S contain the singular values of A arranged in
//    descending order of magnitude.
//
//    Output, double E[MM], where MM = min(M+1,N), ordinarily contains zeros.
//    However see the discussion of INFO for exceptions.
//
//    Output, double U[LDU*K].  If JOBA = 1 then K = M;
//    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of left singular
//    vectors.  U is not referenced if JOBA = 0.  If M <= N or if JOBA = 2, then
//    U may be identified with A in the subroutine call.
//
//    Input, int LDU, the leading dimension of the array U.
//    LDU must be at least M.
//
//    Output, double V[LDV*N], the N by N matrix of right singular vectors.
//    V is not referenced if JOB is 0.  If N <= M, then V may be identified
//    with A in the subroutine call.
//
//    Input, int LDV, the leading dimension of the array V.
//    LDV must be at least N.
//
//    Workspace, double WORK[M].
//
//    Input, int JOB, controls the computation of the singular
//    vectors.  It has the decimal expansion AB with the following meaning:
//      A =  0, do not compute the left singular vectors.
//      A =  1, return the M left singular vectors in U.
//      A >= 2, return the first min(M,N) singular vectors in U.
//      B =  0, do not compute the right singular vectors.
//      B =  1, return the right singular vectors in V.
//
//    Output, int *DSVDC, status indicator INFO.
//    The singular values (and their corresponding singular vectors)
//    S(*INFO+1), S(*INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
//    Thus if *INFO is 0, all the singular values and their vectors are
//    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
//    matrix with the elements of S on its diagonal and the elements of E on
//    its superdiagonal.  Thus the singular values of A and B are the same.
//
{
  double b;
  double c;
  double cs;
  double el;
  double emm1;
  double f;
  double g;
  int i;
  int info;
  int iter;
  int j;
  int jobu;
  int k;
  int kase;
  int kk;
  int l;
  int ll;
  int lls;
  int ls;
  int lu;
  int maxit = 30;
  int mm;
  int mm1;
  int mn;
  int nct;
  int nctp1;
  int ncu;
  int nrt;
  int nrtp1;
  double scale;
  double shift;
  double sl;
  double sm;
  double smm1;
  double sn;
  double t;
  double t1;
  double test;
  bool wantu;
  bool wantv;
  double ztest;
//
//  Determine what is to be computed.
//
  info = 0;
  wantu = false;
  wantv = false;
  jobu = ( job % 100 ) / 10;

  if ( 1 < jobu )
  {
    ncu = i4_min ( m, n );
  }
  else
  {
    ncu = m;
  }

  if ( jobu != 0 )
  {
    wantu = true;
  }

  if ( ( job % 10 ) != 0 )
  {
    wantv = true;
  }
//
//  Reduce A to bidiagonal form, storing the diagonal elements
//  in S and the super-diagonal elements in E.
//
  nct = i4_min ( m-1, n );
  nrt = i4_max ( 0, i4_min ( m, n-2 ) );
  lu = i4_max ( nct, nrt );

  for ( l = 1; l <= lu; l++ )
  {
//
//  Compute the transformation for the L-th column and
//  place the L-th diagonal in S(L).
//
    if ( l <= nct )
    {
      s[l-1] = dnrm2 ( m-l+1, a+l-1+(l-1)*lda, 1 );

      if ( s[l-1] != 0.0 )
      {
        if ( a[l-1+(l-1)*lda] != 0.0 )
        {
          s[l-1] = r8_sign ( a[l-1+(l-1)*lda] ) * fabs ( s[l-1] );
        }
        dscal ( m-l+1, 1.0 / s[l-1], a+l-1+(l-1)*lda, 1 );
        a[l-1+(l-1)*lda] = 1.0 + a[l-1+(l-1)*lda];
      }
      s[l-1] = -s[l-1];
    }

    for ( j = l+1; j <= n; j++ )
    {
//
//  Apply the transformation.
//
      if ( l <= nct && s[l-1] != 0.0 )
      {
        t = - ddot ( m-l+1, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 ) 
          / a[l-1+(l-1)*lda];
        daxpy ( m-l+1, t, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 );
      }
//
//  Place the L-th row of A into E for the
//  subsequent calculation of the row transformation.
//
      e[j-1] = a[l-1+(j-1)*lda];
    }
//
//  Place the transformation in U for subsequent back multiplication.
//
    if ( wantu && l <= nct )
    {
      for ( i = l; i <= m; i++ )
      {
        u[i-1+(l-1)*ldu] = a[i-1+(l-1)*lda];
      }
    }

    if ( l <= nrt )
    {
//
//  Compute the L-th row transformation and place the
//  L-th superdiagonal in E(L).
//
      e[l-1] = dnrm2 ( n-l, e+l, 1 );

      if ( e[l-1] != 0.0 )
      {
        if ( e[l] != 0.0 )
        {
          e[l-1] = r8_sign ( e[l] ) * fabs ( e[l-1] );
        }
        dscal ( n-l, 1.0 / e[l-1], e+l, 1 );
        e[l] = 1.0 + e[l];
      }

      e[l-1] = -e[l-1];
//
//  Apply the transformation.
//
      if ( l+1 <= m && e[l-1] != 0.0 )
      {
        for ( j = l+1; j <= m; j++ )
        {
          work[j-1] = 0.0;
        }

        for ( j = l+1; j <= n; j++ )
        {
          daxpy ( m-l, e[j-1], a+l+(j-1)*lda, 1, work+l, 1 );
        }

        for ( j = l+1; j <= n; j++ )
        {
          daxpy ( m-l, -e[j-1]/e[l], work+l, 1, a+l+(j-1)*lda, 1 );
        }
      }
//
//  Place the transformation in V for subsequent back multiplication.
//
      if ( wantv )
      {
        for ( j = l+1; j <= n; j++ )
        {
          v[j-1+(l-1)*ldv] = e[j-1];
        }
      }
    }
  }
//
//  Set up the final bidiagonal matrix of order MN.
//
  mn = i4_min ( m + 1, n );
  nctp1 = nct + 1;
  nrtp1 = nrt + 1;

  if ( nct < n )
  {
    s[nctp1-1] = a[nctp1-1+(nctp1-1)*lda];
  }

  if ( m < mn )
  {
    s[mn-1] = 0.0;
  }

  if ( nrtp1 < mn )
  {
    e[nrtp1-1] = a[nrtp1-1+(mn-1)*lda];
  }

  e[mn-1] = 0.0;
//
//  If required, generate U.
//
  if ( wantu )
  {
    for ( i = 1; i <= m; i++ )
    {
      for ( j = nctp1; j <= ncu; j++ )
      {
        u[(i-1)+(j-1)*ldu] = 0.0;
      }
    }

    for ( j = nctp1; j <= ncu; j++ )
    {
      u[j-1+(j-1)*ldu] = 1.0;
    }

    for ( ll = 1; ll <= nct; ll++ )
    {
      l = nct - ll + 1;

      if ( s[l-1] != 0.0 )
      {
        for ( j = l+1; j <= ncu; j++ )
        {
          t = - ddot ( m-l+1, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 ) 
            / u[l-1+(l-1)*ldu];
          daxpy ( m-l+1, t, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 );
        }

        dscal ( m-l+1, -1.0, u+(l-1)+(l-1)*ldu, 1 );
        u[l-1+(l-1)*ldu] = 1.0 + u[l-1+(l-1)*ldu];
        for ( i = 1; i <= l-1; i++ )
        {
          u[i-1+(l-1)*ldu] = 0.0;
        }
      }
      else
      {
        for ( i = 1; i <= m; i++ )
        {
          u[i-1+(l-1)*ldu] = 0.0;
        }
        u[l-1+(l-1)*ldu] = 1.0;
      }
    }
  }
//
//  If it is required, generate V.
//
  if ( wantv )
  {
    for ( ll = 1; ll <= n; ll++ )
    {
      l = n - ll + 1;

      if ( l <= nrt && e[l-1] != 0.0 )
      {
        for ( j = l+1; j <= n; j++ )
        {
          t = - ddot ( n-l, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 ) 
            / v[l+(l-1)*ldv];
          daxpy ( n-l, t, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 );
        }

      }
      for ( i = 1; i <= n; i++ )
      {
        v[i-1+(l-1)*ldv] = 0.0;
      }
      v[l-1+(l-1)*ldv] = 1.0;
    }
  }
//
//  Main iteration loop for the singular values.
//
  mm = mn;
  iter = 0;

  while ( 0 < mn )
  {
//
//  If too many iterations have been performed, set flag and return.
//
    if ( maxit <= iter )
    {
      info = mn;
      return info;
    }
//
//  This section of the program inspects for
//  negligible elements in the S and E arrays.
//
//  On completion the variables KASE and L are set as follows:
//
//  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
//  KASE = 2     if S(L) is negligible and L < MN
//  KASE = 3     if E(L-1) is negligible, L < MN, and
//               S(L), ..., S(MN) are not negligible (QR step).
//  KASE = 4     if E(MN-1) is negligible (convergence).
//
    for ( ll = 1; ll <= mn; ll++ )
    {
      l = mn - ll;

      if ( l == 0 )
      {
        break;
      }

      test = fabs ( s[l-1] ) + fabs ( s[l] );
      ztest = test + fabs ( e[l-1] );

      if ( ztest == test )
      {
        e[l-1] = 0.0;
        break;
      }
    }

    if ( l == mn - 1 )
    {
      kase = 4;
    }
    else
    {
      for ( lls = l + 1; lls <= mn + 1; lls++ )
      {
        ls = mn - lls + l + 1;

        if ( ls == l )
        {
          break;
        }

        test = 0.0;
        if ( ls != mn )
        {
          test = test + fabs ( e[ls-1] );
        }

        if ( ls != l + 1 )
        {
          test = test + fabs ( e[ls-2] );
        }

        ztest = test + fabs ( s[ls-1] );

        if ( ztest == test )
        {
          s[ls-1] = 0.0;
          break;
        }

      }

      if ( ls == l )
      {
        kase = 3;
      }
      else if ( ls == mn )
      {
        kase = 1;
      }
      else
      {
        kase = 2;
        l = ls;
      }
    }

    l = l + 1;
//
//  Deflate negligible S(MN).
//
    if ( kase == 1 )
    {
      mm1 = mn - 1;
      f = e[mn-2];
      e[mn-2] = 0.0;

      for ( kk = 1; kk <= mm1; kk++ )
      {
        k = mm1 - kk + l;
        t1 = s[k-1];
        drotg ( &t1, &f, &cs, &sn );
        s[k-1] = t1;

        if ( k != l )
        {
          f = -sn * e[k-2];
          e[k-2] = cs * e[k-2];
        }

        if ( wantv )
        {
          drot ( n, v+0+(k-1)*ldv, 1, v+0+(mn-1)*ldv, 1, cs, sn );
        }
      }
    }
//
//  Split at negligible S(L).
//
    else if ( kase == 2 )
    {
      f = e[l-2];
      e[l-2] = 0.0;

      for ( k = l; k <= mn; k++ )
      {
        t1 = s[k-1];
        drotg ( &t1, &f, &cs, &sn );
        s[k-1] = t1;
        f = - sn * e[k-1];
        e[k-1] = cs * e[k-1];
        if ( wantu )
        {
          drot ( m, u+0+(k-1)*ldu, 1, u+0+(l-2)*ldu, 1, cs, sn );
        }
      }
    }
//
//  Perform one QR step.
//
    else if ( kase == 3 )
    {
//
//  Calculate the shift.
//
      scale = r8_max ( fabs ( s[mn-1] ), 
              r8_max ( fabs ( s[mn-2] ), 
              r8_max ( fabs ( e[mn-2] ), 
              r8_max ( fabs ( s[l-1] ), fabs ( e[l-1] ) ) ) ) );

      sm = s[mn-1] / scale;
      smm1 = s[mn-2] / scale;
      emm1 = e[mn-2] / scale;
      sl = s[l-1] / scale;
      el = e[l-1] / scale;
      b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0;
      c = ( sm * emm1 ) * ( sm * emm1 );
      shift = 0.0;

      if ( b != 0.0 || c != 0.0 )
      {
        shift = sqrt ( b * b + c );
        if ( b < 0.0 )
        {
          shift = -shift;
        }
        shift = c / ( b + shift );
      }

      f = ( sl + sm ) * ( sl - sm ) - shift;
      g = sl * el;
//
//  Chase zeros.
//
      mm1 = mn - 1;

      for ( k = l; k <= mm1; k++ )
      {
        drotg ( &f, &g, &cs, &sn );

        if ( k != l )
        {
          e[k-2] = f;
        }

        f = cs * s[k-1] + sn * e[k-1];
        e[k-1] = cs * e[k-1] - sn * s[k-1];
        g = sn * s[k];
        s[k] = cs * s[k];

        if ( wantv )
        {
          drot ( n, v+0+(k-1)*ldv, 1, v+0+k*ldv, 1, cs, sn );
        }

        drotg ( &f, &g, &cs, &sn );
        s[k-1] = f;
        f = cs * e[k-1] + sn * s[k];
        s[k] = -sn * e[k-1] + cs * s[k];
        g = sn * e[k];
        e[k] = cs * e[k];

        if ( wantu && k < m )
        {
          drot ( m, u+0+(k-1)*ldu, 1, u+0+k*ldu, 1, cs, sn );
        }
      }
      e[mn-2] = f;
      iter = iter + 1;
    }
//
//  Convergence.
//
    else if ( kase == 4 )
    {
//
//  Make the singular value nonnegative.
//
      if ( s[l-1] < 0.0 )
      {
        s[l-1] = -s[l-1];
        if ( wantv )
        {
          dscal ( n, -1.0, v+0+(l-1)*ldv, 1 );
        }
      }
//
//  Order the singular value.
//
      for ( ; ; )
      {
        if ( l == mm )
        {
          break;
        }

        if ( s[l] <= s[l-1] )
        {
          break;
        }

        t = s[l-1];
        s[l-1] = s[l];
        s[l] = t;

        if ( wantv && l < n )
        {
          dswap ( n, v+0+(l-1)*ldv, 1, v+0+l*ldv, 1 );
        }

        if ( wantu && l < m )
        {
          dswap ( m, u+0+(l-1)*ldu, 1, u+0+l*ldu, 1 );
        }

        l = l + 1;
      }
      iter = 0;
      mn = mn - 1;
    }
  }

  return info;
}
//****************************************************************************80

void dswap ( int n, double x[], int incx, double y[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DSWAP interchanges two vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double X[*], one of the vectors to swap.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], one of the vectors to swap.
//
//    Input, int INCY, the increment between successive elements of Y.
//
{
  int i;
  int ix;
  int iy;
  int m;
  double temp;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 && incy == 1 )
  {
    m = n % 3;

    for ( i = 0; i < m; i++ )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;
    }

    for ( i = m; i < n; i = i + 3 )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;

      temp = x[i+1];
      x[i+1] = y[i+1];
      y[i+1] = temp;

      temp = x[i+2];
      x[i+2] = y[i+2];
      y[i+2] = temp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      temp = x[ix];
      x[ix] = y[iy];
      y[iy] = temp;
      ix = ix + incx;
      iy = iy + incy;
    }

  }

  return;
}
//****************************************************************************80

double *normal_solve ( int m, int n, double a[], double b[], int &flag )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_SOLVE solves a linear system using the normal equations.
//
//  Discussion:
//
//    Given a presumably rectangular MxN system of equations A*x=b, this routine
//    sets up the NxN system A'*A*x=A'b.  Assuming N <= M, and that A has full
//    column rank, the system will be solvable, and the vector x that is returned
//    will minimize the Euclidean norm of the residual.
//
//    One drawback to this approach is that the condition number of the linear
//    system A'*A is effectively the square of the condition number of A, 
//    meaning that there is a substantial loss of accuracy.
//
//    Thanks to David Doria for pointing out that this procedure was missing
//    "delete[]" statements at the end, 08 April 2013.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//    It must be the case that N <= M.
//
//    Input, double A[M*N], the matrix.
//    The matrix must have full column rank.
//
//    Input, double B[M], the right hand side.
//
//    Output, int &FLAG,
//    0, no error was detected.
//    1, an error occurred.
//
//    Output, double NORMAL_SOLVE[N], the least squares solution.
//
{
  double *at;
  double *ata;
  double *ata_c;
  double *atb;
  double *x;

  flag = 0;

  if ( m < n )
  {
    flag = 1;
    return NULL;
  }

  at = r8mat_transpose_new ( m, n, a );

  ata = r8mat_mm_new ( n, m, n, at, a );

  ata_c = r8mat_cholesky_factor ( n, ata, flag );

  if ( flag != 0 )
  {
    return NULL;
  }

  atb = r8mat_mv_new ( n, m, at, b );

  x = r8mat_cholesky_solve ( n, ata_c, atb );

  delete [] at;
  delete [] ata;
  delete [] ata_c;
  delete [] atb;

  return x;
}
//****************************************************************************80

double *qr_solve ( int m, int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    QR_SOLVE solves a linear system in the least squares sense.
//
//  Discussion:
//
//    If the matrix A has full column rank, then the solution X should be the
//    unique vector that minimizes the Euclidean norm of the residual.
//
//    If the matrix A does not have full column rank, then the solution is
//    not unique; the vector X will minimize the residual norm, but so will
//    various other vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input, double A[M*N], the matrix.
//
//    Input, double B[M], the right hand side.
//
//    Output, double QR_SOLVE[N], the least squares solution.
//
{
  double *a_qr;
  int itask;
  int *jpvt;
  int kr;
  int lda;
  double *qraux;
  double *r;
  double tol;
  double *x;

  a_qr = r8mat_copy_new ( m, n, a );
  lda = m;
  tol = r8_epsilon ( ) / r8mat_amax ( m, n, a_qr );
  x = new double[n];
  jpvt = new int[n];
  qraux = new double[n];
  r = new double[m];
  itask = 1;

  dqrls ( a_qr, lda, m, n, tol, kr, b, x, r, jpvt, qraux, itask );

  delete [] a_qr;
  delete [] jpvt;
  delete [] qraux;
  delete [] r;

  return x;
}
//****************************************************************************80

double *svd_solve ( int m, int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    SVD_SOLVE solves a linear system in the least squares sense.
//
//  Discussion:
//
//    The vector X returned by this routine should always minimize the 
//    Euclidean norm of the residual ||A*x-b||.
//
//    If the matrix A does not have full column rank, then there are multiple
//    vectors that attain the minimum residual.  In that case, the vector
//    X returned by this routine is the unique such minimizer that has the 
//    the minimum possible Euclidean norm, that is, ||A*x-b|| and ||x||
//    are both minimized.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input, double A[M*N], the matrix.
//
//    Input, double B[M], the right hand side.
//
//    Output, double SVD_SOLVE[N], the least squares solution.
//
{
  double *a_copy;
  double *e;
  int i;
  int info;
  int lda;
  int ldu;
  int ldv;
  int job;
  double *sdiag;
  double smax;
  double stol;
  double *sub;
  double *u;
  double *ub;
  double *v;
  double *work;
  double *x;
//
//  Get the SVD.
//
  a_copy = r8mat_copy_new ( m, n, a );
  lda = m;
  sdiag = new double[i4_max ( m + 1, n )];
  e = new double[i4_max ( m + 1, n )];
  u = new double[m*m];
  ldu = m;
  v = new double[n*n];
  ldv = n;
  work = new double[m];
  job = 11;

  info = dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job );

  delete [] work;

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "SVD_SOLVE - Failure!\n";
    cerr << "  The SVD could not be calculated.\n";
    cerr << "  LINPACK routine DSVDC returned a nonzero\n";
    cerr << "  value of the error flag, INFO = " << info << "\n";
    exit ( 1 );
  }

  ub = r8mat_mtv_new ( m, m, u, b );
//
//  For singular problems, there may be tiny but nonzero singular values
//  that should be ignored.  This is a reasonable attempt to avoid such 
//  problems, although in general, the user might wish to control the tolerance.
//
  smax = r8vec_max ( n, sdiag );
  if ( smax <= r8_epsilon ( ) )
  {
    smax = 1.0;
  }

  stol = r8_epsilon ( ) * smax;

  sub = new double[n];

  for ( i = 0; i < n; i++ )
  {
    sub[i] = 0.0;
    if ( i < m )
    {
      if ( stol <= sdiag[i] )
      {
        sub[i] = ub[i] / sdiag[i];
      }
    }
  }

  x = r8mat_mv_new ( n, n, v, sub );

  delete [] a_copy;
  delete [] e;
  delete [] sdiag;
  delete [] sub;
  delete [] u;
  delete [] ub;
  delete [] v;

  return x;
}
//Wrapper function by Nick
//Solves a matrix system Mx=y for x, where 
//M is nRow x Ncol, 
//y is nCol x nSol, and 
//x is nRow x nSol
void qr_solve_mat(int nRow, int nCol, int nSol, double* M, double* y, double* x) {
    int i;
    double* qrResult;
    for (i = 0; i < nSol; i++) {
        qrResult = qr_solve(nRow, nCol, M, &y[i * nRow]);
        memcpy_s(&x[i * nCol], nCol * sizeof(double), qrResult, nCol * sizeof(double));
        free(qrResult);
    }
}