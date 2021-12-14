#include "utils.h"

double *rand_dvec(size_t m, double lower, double upper) {
  double *res = malloc(m * sizeof(double));
  for (size_t i = 0; i < m; i++) {
    res[i] = lower + (upper - lower) * ((double) rand() / RAND_MAX);
  }
  return res;
}

void print_dvec(double* vec, size_t m) {
  for (size_t i = 0; i < m; i++) {
    printf("%lf ", vec[i]);
  }
  printf("\n");
}

void print_svec(float* vec, size_t m) {
  for (size_t i = 0; i < m; i++) {
    printf("%f ", vec[i]);
  }
  printf("\n");
}

void print_dmat(double* mat, size_t n, size_t m) {
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      printf("%lf ", mat[i * m + j]);
    }
    printf("\n");
  }
}