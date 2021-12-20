#include "poisson1D.h"
#include "blaslapack_headers.h"
#include "utils.h"

void drich(double omega, double *A, double *x, double *b, size_t n, size_t m) {

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

  drich(0.5, A, x, b, n, m);
  free(A);
  return 0;
}