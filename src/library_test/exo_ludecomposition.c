
#include "poisson1D.h"
#include "blaslapack_headers.h"
#include "utils.h"
#include <time.h>

// LU decomposition of a tridiagonal matrix in general band format
// Double precision version
void dtridb_lude(double *A, size_t la) {
  // Aliases
  double *c = A;
  double *b = A + la;
  double *a = A + 2 * la;

  // Optimized lu decomposition for tri-diag matrices
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

void run_lude(size_t n) {
  double *A = malloc(n * 3 * sizeof(double));
  makeRowMajorGBand(A, 3, n, 0);

  dtridb_lude(A, n);

  free(A);
}

void run_dgbtrf(size_t n) {
  double *A = malloc(n * 4 * sizeof(double));
  int* ipiv = malloc(n * sizeof(int));
  makeRowMajorGBand(A, 3, n, 1);

  LAPACKE_dgbtrf(LAPACK_ROW_MAJOR, n, n, 1, 1, A, n, ipiv);

  free(ipiv);
  free(A);
}

int main() {
  FILE* out = fopen("lude_time.dat", "w");
  const int kSampleSize = 5;

  for (int i = 10; i <= 5000; i++) {
    clock_t begin = clock();
    // This is a really simple time measurement
    // since we take into account the matrix creation
    // but this is easier
    for (size_t j = 0; j < kSampleSize; j++) {
      run_lude(i);
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC / kSampleSize;
    fprintf(out, "%i %.8f\n", i, time_spent);
  }
  fclose(out);


  out = fopen("dgbtrf_time.dat", "w");

  for (int i = 10; i <= 5000; i++) {
    clock_t begin = clock();
    // This is a really simple time measurement
    // since we take into account the matrix creation
    // but this is easier
    for (size_t j = 0; j < kSampleSize; j++) {
      run_dgbtrf(i);
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC / kSampleSize;
    fprintf(out, "%i %.8f\n", i, time_spent);
  }
  fclose(out);

  return 0;
}