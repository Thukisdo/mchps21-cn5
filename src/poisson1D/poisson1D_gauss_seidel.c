/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "poisson1D.h"
#include <time.h>

void run_gseidel(int nbpoints) {
  int ierr;
  int jj;
  int la;
  int ku, kl, lab, kv;
  // Initialisation ???
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  // Initialisation ???
  double *RHS, *SOL, *EX_SOL, *X;
  // Initialisation ???
  double *AB;

  double temp, relres;

  double opt_alpha;

  /* Size of the problem */
  NRHS = 1;
  la = nbpoints - 2;

  /* Dirichlet Boundary conditions */
  T0 = -5.0;
  T1 = 5.0;

  // Error checking ???
  RHS = (double *) malloc(sizeof(double) * la);
  SOL = (double *) calloc(la, sizeof(double));
  EX_SOL = (double *) malloc(sizeof(double) * la);
  X = (double *) malloc(sizeof(double) * la);

  cblas_dscal(la, 0.0, SOL, 1);

  /* Setup the Poisson 1D problem */
  /* General Band Storage */
  setGridPoints(X, &la);
  set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
  computeAnalyticalSolution(EX_SOL, X, &la, &T0, &T1);

  writeVec(RHS, &la, "RHS.dat");
  writeVec(EX_SOL, &la, "EX_SOL.dat");
  writeVec(X, &la, "X_grid.dat");

  kv = 0;
  ku = 1;
  kl = 1;
  lab = kv + kl + ku + 1;

  // Error checking ???
  AB = (double *) malloc(sizeof(double) * lab * la);
  makeColMajorGBand(AB, lab, la, kv);

  /* Solve */
  double tol = 1e-8;
  int maxit = 10000;
  gaussSeidel(AB, RHS, SOL, lab, la, ku, kl, tol, maxit);

  /* Write solution */
  writeVec(SOL, &la, "SOL.dat");

  free(RHS);
  free(SOL);
  free(EX_SOL);
  free(X);
  free(AB);
  // ??? Return ???
}

int main(int argc, char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  FILE* out = fopen("gs_time.dat", "w");
  const int kSampleSize = 5;

  for (int i = 10; i <= 300; i++) {
    clock_t begin = clock();
    // This is a really simple time measurement
    // since we take into account the matrix creation
    // but this is easier
    for (size_t j = 0; j < kSampleSize; j++) {
      run_gseidel(i);
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC / kSampleSize;
    fprintf(out, "%i %f\n", i, time_spent);
  }
  fclose(out);
}
