
#include "poisson1D.h"
#include "blaslapack_headers.h"
#include "utils.h"

// LU decomposition of a tridiagonal matrix in general band format
// Double precision version
void dtridb_lude(double *A, size_t la) {
  // Aliases
  double *c = A;
  double *b = A + la;
  double *a = A + 2 * la;

  // Optimized tri-diag lu decomposition
  // Since we use a general band storage,
  // We must be careful on how we index our elements
  for (size_t i = 1; i < la; i++) {
    // a0 = a[0] && bi = b[i]
    double lk = a[i - 1] / b[i - 1];
    // c0 = c[0 + 1]
    b[i] -= lk * c[i];
    a[i - 1] = lk;
  }

}

int main() {
  size_t n = 12;
  double *A = malloc(n * 3 * sizeof(double));

  set_GB_operator_rowMajor_poisson1D(A, 3, n, 0);
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < n; j++) {
      printf("%f ", A[i * n + j]);
    }
    printf("\n");
  }
  printf("\n");


  dtridb_lude(A, n);
  print_dmat(A, 3, n);

  free(A);
  return 0;
}