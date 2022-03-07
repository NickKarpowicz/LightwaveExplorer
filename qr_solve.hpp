void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy );
double ddot ( int n, double dx[], int incx, double dy[], int incy );
double dnrm2 ( int n, double x[], int incx );
void dqrank ( double a[], int lda, int m, int n, double tol, int &kr, 
  int jpvt[], double qraux[] );
void dqrdc ( double a[], int lda, int n, int p, double qraux[], int jpvt[], 
  double work[], int job );
int dqrls ( double a[], int lda, int m, int n, double tol, int &kr, double b[], 
  double x[], double rsd[], int jpvt[], double qraux[], int itask );
void dqrlss ( double a[], int lda, int m, int n, int kr, double b[], double x[], 
  double rsd[], int jpvt[], double qraux[] );
int dqrsl ( double a[], int lda, int n, int k, double qraux[], double y[], 
  double qy[], double qty[], double b[], double rsd[], double ab[], int job );
void drot ( int n, double x[], int incx, double y[], int incy, double c,
  double s );
void drotg ( double *sa, double *sb, double *c, double *s );
void dscal ( int n, double sa, double x[], int incx );
int dsvdc ( double a[], int lda, int m, int n, double s[], double e[], 
  double u[], int ldu, double v[], int ldv, double work[], int job );
void dswap ( int n, double x[], int incx, double y[], int incy );
double *normal_solve ( int m, int n, double a[], double b[], int &flag );
double *qr_solve ( int m, int n, double a[], double b[] );
double *svd_solve ( int m, int n, double a[], double b[] );
void qr_solve_mat(int nRow, int nCol, int nSol, double* M, double* y, double* x);