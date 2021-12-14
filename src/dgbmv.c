#include "blaslapack_headers.h"
#include "lib_poisson1D.h"
#include "clapack.h"
#include "utils.h"

// Performs a matrix-vector product using the DGBMV BLAS routine.
// y = alpha * A * x + beta * y
void dgbmv_row_major(double alpha, double *A, double *x, int la, int lab, double beta, double *y) {

  cblas_dgbmv(/*Order*/CblasRowMajor,/*Trans*/ CblasNoTrans,/*A rows*/ la, /*A cols*/
                       la, /*N lower diagonals*/ 1, /*N Upper diagonals*/ 1, alpha, A, lab, x, /*inc x*/ 1,
                       beta, y, /*inc y*/1);
}

// Performs a matrix-vector product using the DGBMV BLAS routine.
// y = alpha * A * x + beta * y
void dgbmv_col_major(double alpha, double *A, double *x, int la, int lab, double beta, double *y) {
  cblas_dgbmv(/*Order*/CblasColMajor,/*Trans*/ CblasNoTrans,/*A rows*/ la, /*A cols*/
                       la, /*N lower diagonals*/ 1, /*N Upper diagonals*/ 1, alpha, A, lab, x, /*inc x*/ 1,
                       beta, y, /*inc y*/1);
}

int main() {

  printf("--------- DGBMV 1D ---------\n\n");
  // la = band length
  int nbpoints = 12, la = nbpoints - 2;
  double *A = NULL, *x = NULL, *y = NULL;

  // Number of bands
  int lab = 3;
  A = (double *) malloc(la * lab * sizeof(double));

  x = random_dvec(la, 0, 0);
  y = random_dvec(la, 0, 0);


  // Random initialization
  set_GB_operator_colMajor_poisson1D(A, lab, la, 0);
  for (size_t i = 0; i < lab; i++) {
    for (size_t j = 0; j < la; j++) {
      printf("%f ", A[i * la + j]);
    }
    printf("\n");
  }
  for (size_t i = 0; i < la; i++) {
    x = random_dvec(la, 0, 0);
    x[i] = 1.0;
    printf("i = %zu\n", i);
    dgbmv_col_major(1.0, A, x, la, lab, 1.0, y);
    for (size_t j = 0; j < la; j++) {
      printf("y[%d] = %f\n", j, y[j]);
    }
    printf("\n");
    free(x);
  }
  free(y);

  // Random initialization
  // set_GB_operator_rowMajor_poisson1D(A, lab, la, 0);
  // dgbmv_row_major(1, A, x, la, lab, 1, y);

  return 0;
}