#include "poisson1D.h"
#include "blaslapack_headers.h"
#include "utils.h"

void drich(double omega, double *A, double *x, double *b, size_t n, size_t m) {
  double *work = malloc(sizeof(double) * n);

  double *inv_diag = malloc(sizeof(double) * n);
  for (size_t i = 0; i < n; i++) {
    inv_diag[i] = 1.0 / A[i * n + i];
  }

  double res_norm = 1;
  size_t counter = 0;
  while(res_norm > 1e-2) {
    cblas_dcopy(n, b, 1, work, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, m, -1.0, A, m, x, 1, 1.0, work, 1);
    cblas_daxpy(n, omega, work, 1, x, 1);
    res_norm = cblas_dnrm2(n, work, 1);
    printf("%zu %f\n", counter, res_norm);
    counter++;
  }

  /*printf("Sol:\n");
  for (size_t i = 0; i < n; i++) {
    printf("%.16f ", x[i]);
  }
  printf("\n\n");

  printf("Residu:\n");
  for (size_t i = 0; i < n; i++) {
    printf("%.16f ", work[i]);
  }*/
  free(work);
}


int main() {
  size_t n = 100, m = 100;
  double *A = malloc(sizeof(double) * n * m);
  double *x = malloc(sizeof(double) * n);
  double *b = malloc(sizeof(double) * n);

  for (size_t i = 0; i < n; i++) {
    A[i * m + i] = 2.0;
  }
  for (size_t i = 0; i < n - 1; i++) {
    A[i * m + i + 1] = -1.0;
    A[(i + 1) * m + i] = -1.0;
  }

  for (size_t i = 0; i < n; i++) {
    x[i] = 0.0;
    b[i] = 1.0;
  }

  /*print_dmat(A, n, m);
  printf("\n");
  print_dvec(b, n);
  printf("\n");
  print_dvec(x, n);
  printf("\n");*/


  drich(0.5, A, x, b, n, m);
  free(A);
  return 0;
}