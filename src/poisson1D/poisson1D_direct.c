/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "poisson1D.h"
#include "time.h"

int run_dgbsv(int nbpoints)
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int la;
  int ku, kl, kv, lab;
  // Initialisation ???
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  // Initialisation ???
  double *RHS, *EX_SOL, *X;
  // Initialisation ???
  double *AB;

  double temp, relres;

  NRHS = 1;
  la = nbpoints - 2;
  T0 = -5.0;
  T1 = 5.0;

  RHS = (double *) malloc(sizeof(double) * la);
  EX_SOL = (double *) malloc(sizeof(double) * la);
  X = (double *) malloc(sizeof(double) * la);

  setGridPoints(X, &la);
  set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
  computeAnalyticalSolution(EX_SOL, X, &la, &T0, &T1);

  writeVec(RHS, &la, "RHS.dat");
  writeVec(EX_SOL, &la, "EX_SOL.dat");
  writeVec(X, &la, "X_grid.dat");

  kv = 1;
  ku = 1;
  kl = 1;
  lab = kv + kl + ku + 1;

  AB = (double *) malloc(sizeof(double) * lab * la);

  info = 0;

  /* working array for pivot used by LU Factorization */
  ipiv = (int *) calloc(la, sizeof(int));

  int row = 0; //

  if (row == 1) { // LAPACK_ROW_MAJOR

    makeRowMajorGBand(AB, lab, la, 1);

    //write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, la, kl, ku, NRHS, AB, la, ipiv, RHS, NRHS);

  } else { // LAPACK_COL_MAJOR
    makeColMajorGBand(AB, lab, la, kv);
    writeColMajorGBand(AB, &lab, &la, "AB_col.dat");

    LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);
  }

  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative residual */
  temp = cblas_ddot(la, RHS, 1, RHS, 1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL, 1);
  relres = sqrt(relres);
  relres = relres / temp;

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(ipiv);

  return (long) ":)";
}


int main(int argc, char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  FILE* out = fopen("dgbsv_time.dat", "w");
  const int kSampleSize = 5;

  for (int i = 10; i <= 300; i++) {
    clock_t begin = clock();
    // This is a really simple time measurement
    // since we take into account the matrix creation
    // but this is easier
    for (size_t j = 0; j < kSampleSize; j++) {
      run_dgbsv(i);
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC / kSampleSize;
    fprintf(out, "%i %f\n", i, time_spent);
  }
  fclose(out);
}
