/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "poisson1D.h"
#include <stdio.h>
#include <stdlib.h>


Poisson1DProblem createPoisson1D(int la, double T0, double T1,
                                 int write_to_file) {
  Poisson1DProblem res = {NULL, NULL, NULL, NULL, 0};

  res.rhs = malloc(sizeof(double) * la);
  res.sol = calloc(la, sizeof(double));
  res.exact_sol = malloc(sizeof(double) * la);
  res.x = malloc(sizeof(double) * la);
  res.la = la;

  cblas_dscal(la, 0.0, res.sol, 1);

  setGridPoints(res.x, &la);
  set_dense_RHS_DBC_1D(res.rhs, &la, &T0, &T1);
  computeAnalyticalSolution(res.exact_sol, res.x, &la, &T0, &T1);

  if (write_to_file) {
    writeVec(res.rhs, &la, "RHS.dat");
    writeVec(res.exact_sol, &la, "EX_SOL.dat");
    writeVec(res.x, &la, "X_grid.dat");
  }

  return res;
}

void freePoisson1D(Poisson1DProblem p) {}

void makeRowMajorGBand(double *AB, int lab, int la, int kv) {

  for (int j = 0; j < la * kv; j++) { AB[j] = 0; }

  for (int j = 0; j < la; j++) { AB[j + (kv) *la] = -1.0; }

  for (int j = 0; j < la; j++) { AB[j + (kv + 1) * la] = 2.0; }

  for (int j = 0; j < la; j++) { AB[j + (kv + 2) * la] = -1.0; }

  AB[0] = 0;
  AB[lab * la - 1] = 0;
}

void makeColMajorGBand(double *AB, int lab, int la, int kv) {
  for (int j = 0; j < la; j++) {
    int k = j * lab;
    if (kv >= 0) {
      for (int i = 0; i < kv; i++) { AB[k + i] = 0.0; }
    }
    AB[k + kv] = -1.0;
    AB[k + kv + 1] = 2.0;
    AB[k + kv + 2] = -1.0;
  }
  AB[0] = 0.0;
  if (kv == 1) { AB[1] = 0; }

  AB[lab * la - 1] = 0.0;
}

void makeColMajorGBIdentity(double *AB, int lab, int la, int kv) {
  for (int j = 0; j < la; j++) {
    int k = j * lab;
    if (kv >= 0) {
      for (int i = 0; i < kv; i++) { AB[k + i] = 0.0; }
    }
    AB[k + kv] = 0.0;
    AB[k + kv + 1] = 1.0;
    AB[k + kv + 2] = 0.0;
  }
  AB[1] = 0.0;
  AB[lab * la - 1] = 0.0;
}

void set_dense_RHS_DBC_1D(double *RHS, int *la, double *BC0, double *BC1) {
  int jj;
  RHS[0] = *BC0;
  RHS[(*la) - 1] = *BC1;
  for (jj = 1; jj < (*la) - 1; jj++) { RHS[jj] = 0.0; }
}

void computeAnalyticalSolution(double *EX_SOL, double *X, int *la, double *BC0,
                               double *BC1) {
  int jj;
  double h, DELTA_T;
  DELTA_T = (*BC1) - (*BC0);
  for (jj = 0; jj < (*la); jj++) { EX_SOL[jj] = (*BC0) + X[jj] * DELTA_T; }
}

void setGridPoints(double *x, int *la) {
  int jj;
  double h;
  h = 1.0 / (1.0 * ((*la) + 1));
  for (jj = 0; jj < (*la); jj++) { x[jj] = (jj + 1) * h; }
}

void writeRowMajorGBand(double *AB, int *lab, int *la, char *filename) {
  FILE *file;
  int ii, jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL) {
    for (ii = 0; ii < (*lab); ii++) {
      for (jj = 0; jj < (*la); jj++) {
        fprintf(file, "%lf\t", AB[ii * (*la) + jj]);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  } else {
    perror(filename);
  }
}

void writeColMajorGBand(double *AB, int *lab, int *la, char *filename) {
  //TODO
}

void writeVec(double *vec, int *la, char *filename) {
  int jj;
  FILE *file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL) {
    for (jj = 0; jj < (*la); jj++) { fprintf(file, "%lf\n", vec[jj]); }
    fclose(file);
  } else {
    perror(filename);
  }
}

void write_xy(double *vec, double *x, int *la, char *filename) {
  int jj;
  FILE *file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL) {
    for (jj = 0; jj < (*la); jj++) {
      fprintf(file, "%lf\t%lf\n", x[jj], vec[jj]);
    }
    fclose(file);
  } else {
    perror(filename);
  }
}

void computeEigenValues(double *eigval, int *la) {
  int ii;
  double scal;
  for (ii = 0; ii < *la; ii++) {
    scal = (1.0 * ii + 1.0) * M_PI_2 * (1.0 / (*la + 1));
    eigval[ii] = sin(scal);
    eigval[ii] = 4 * eigval[ii] * eigval[ii];
  }
}

double computeMaxEigenValue(int la) {
  double eigmax = sin(la * M_PI_2 * (1.0 / (la + 1)));
  eigmax = 4.0 * eigmax * eigmax;
  return eigmax;
}

double computeMinEigenValue(int la) {
  double eigmin = sin(M_PI_2 * (1.0 / (la + 1)));
  eigmin = 4.0 * eigmin * eigmin;
  return eigmin;
}

double computeRichardsonOptAlpha(int la) {
  return 2.0 / (computeMaxEigenValue(la) + computeMinEigenValue(la));
}

void gaussSeidel(double *AB, Poisson1DProblem p, int lab, int ku, int kl,
                 double tol, int maxit, FILE *out) {
  double *work = malloc(sizeof(double) * p.la);

  double res_norm = 1;
  size_t counter = 0;

  // This is suboptimal since M is a lower triangular
  // But I am too lazy to do this properly
  // The inverse of a tridiagonal matrix isn't tridiagonal
  // So we'll use a dense matrix to store the result (a lower triangular)
  double *M = malloc(sizeof(double) * p.la * p.la);

  for (size_t i = 0; i < p.la * p.la; i++) { M[i] = 0.0; }

  // Build the lower diagonal matrix
  for (size_t i = 0; i < p.la; i++) { M[i * p.la + i] = AB[lab * i + 1]; }

  for (size_t i = 0; i < p.la - 1; i++) {
    M[(i + 1) * p.la + i] = AB[lab * i + 2];
  }

  // Inverse it
  int *ipiv = malloc(sizeof(int) * p.la);
  LAPACKE_dgetrf(CblasRowMajor, p.la, p.la, M, p.la, ipiv);
  LAPACKE_dgetri(CblasRowMajor, p.la, M, p.la, ipiv);
  free(ipiv);

  // x(k+1) = x(k) + M^(-1)(b - Ax(k))
  while (res_norm > tol && counter < maxit) {

    cblas_dcopy(p.la, p.rhs, 1, work, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, p.la, p.la, kl, ku, -1.0, AB, lab,
                p.x, 1, 1.0, work, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, p.la, p.la, 1.0, M, p.la, work, 1,
                1.0, p.x, 1);
    if (out) {
      res_norm = cblas_dnrm2(p.la, work, 1);
      res_norm /= cblas_dnrm2(p.la, p.x, 1);
      fprintf(out, "%zu %.16lf\n", counter, res_norm);
    }
    counter++;
  }
  free(M);
  free(work);
}

void richardsonWithAlpha(double *AB, Poisson1DProblem p, int lab, int ku,
                         int kl, double tol, int maxit, FILE *out) {
  double *work = malloc(sizeof(double) * p.la);

  double res_norm = 1;
  size_t counter = 0;

  double alpha_rich = computeRichardsonOptAlpha(p.la);

  // x(k+1) = x(k) + alpha(b - Ax(k))
  while (res_norm > tol && counter < maxit) {

    cblas_dcopy(p.la, p.rhs, 1, work, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, p.la, p.la, kl, ku, -1.0, AB, lab, p.x,
                1, 1.0, work, 1);
    cblas_daxpy(p.la, alpha_rich, work, 1, p.x, 1);
    if (out) {
      res_norm = cblas_dnrm2(p.la, work, 1);
      res_norm /= cblas_dnrm2(p.la, p.x, 1);
      fprintf(out, "%zu %.16lf\n", counter, res_norm);
    }

    counter++;
  }
  free(work);
}
