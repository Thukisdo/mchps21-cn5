/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "poisson1D.h"
#include <time.h>

double runGSeidel(int nbpoints, double tol, int maxiter, FILE* out) {
  int la = nbpoints;
  double T0 = -5.0, T1 = 5;

  int kv = 0;
  int ku = 1;
  int kl = 1;
  int lab = kv + kl + ku + 1;

  Poisson1DProblem p = createPoisson1D(la, T0, T1, 0);

  // Error checking ???
  double *AB = (double *) malloc(sizeof(double) * lab * la);
  makeColMajorGBand(AB, lab, la, kv);

  /* Solve */
  clock_t begin = clock();
  gaussSeidel(AB, p, lab, ku, kl, tol, maxiter, out);
  clock_t end = clock();

  /* Write solution */
  writeVec(p.sol, &la, "SOL.dat");

  free(AB);
  freePoisson1D(p);
  return (double) (end - begin) / CLOCKS_PER_SEC;
}

double perfGSeidel(int nbpoints, double tol, int maxiter) {
  return runGSeidel(nbpoints, tol, maxiter, NULL);
}

void cvGSeidel(int nbpoints, double tol, int maxiter) {
  FILE* f = fopen("gseidel.dat", "w");
  runGSeidel(nbpoints, tol, maxiter, f);
  fclose(f);
}

int main(int argc, char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  FILE *out = fopen("perfGSeidel.dat", "w");
  const int kSampleSize = 5;
  const int maxiter = 10000;

  for (int i = 10; i <= 300; i++) {

    // This is a really simple time measurement
    // since we take into account the matrix creation
    // but this is easier
    double acc = 0;
    for (size_t j = 0; j < kSampleSize; j++) {
      acc += perfGSeidel(i, 10e-8, maxiter);
    }
    acc /= kSampleSize;

    fprintf(out, "%i %f\n", i, acc);
  }
  fclose(out);

  cvGSeidel(50, 10e-8, maxiter);
}
