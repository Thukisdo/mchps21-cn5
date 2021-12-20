/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "poisson1D.h"

void makeRowMajorGBand(double *AB, int lab, int la, int kv) {

  for (int j = 0; j < la * kv; j++) {
    AB[j] = 0;
  }

  for (int j = 0; j < la; j++) {
    AB[j + (kv) * la] = -1.0;
  }

  for (int j = 0; j < la; j++) {
    AB[j + (kv + 1) * la] = 2.0;
  }

  for (int j = 0; j < la; j++) {
    AB[j + (kv + 2) * la] = -1.0;
  }

  AB[0] = 0;
  AB[lab * la - 1] = 0;
}

void makeColMajorGBand(double *AB, int lab, int la, int kv) {
  for (int j = 0; j < la; j++) {
    int k = j * lab;
    if (kv >= 0) {
      for (int i = 0; i < kv; i++) {
        AB[k + i] = 0.0;
      }
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
      for (int i = 0; i < kv; i++) {
        AB[k + i] = 0.0;
      }
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
  for (jj = 1; jj < (*la) - 1; jj++) {
    RHS[jj] = 0.0;
  }
}

void computeAnalyticalSolution(double *EX_SOL, double *X, int *la, double *BC0, double *BC1) {
  int jj;
  double h, DELTA_T;
  DELTA_T = (*BC1) - (*BC0);
  for (jj = 0; jj < (*la); jj++) {
    EX_SOL[jj] = (*BC0) + X[jj] * DELTA_T;
  }
}

void setGridPoints(double *x, int *la) {
  int jj;
  double h;
  h = 1.0 / (1.0 * ((*la) + 1));
  for (jj = 0; jj < (*la); jj++) {
    x[jj] = (jj + 1) * h;
  }
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
    for (jj = 0; jj < (*la); jj++) {
      fprintf(file, "%lf\n", vec[jj]);
    }
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

void gaussSeidel(double *AB, double *RHS, double *X, int lab, int la, int ku, int kl,
                 double tol, int maxit) {
  double *work = malloc(sizeof(double) * la);

  double res_norm = 1;
  size_t counter = 0;

  // This is suboptimal since M is a lower triangular
  // But I am too lazy to do this properly
  // The inverse of a tridiagonal matrix isn't tridiagonal
  // So we'll use a dense matrix to store the result (a lower triangular)
  double *M = malloc(sizeof(double) * la * la);

  for (size_t i = 0; i < la * la; i++) {
    M[i] = 0.0;
  }

  // Build the lower diagonal matrix
  for (size_t i = 0; i < la; i++) {
    M[i * la + i] = AB[lab * i + 1];
  }

  for (size_t i = 0; i < la - 1; i++) {
    M[(i + 1) * la + i] = AB[lab * i + 2];
  }

  // Inverse it
  int* ipiv = malloc(sizeof(int) * la);
  LAPACKE_dgetrf(CblasRowMajor, la, la, M, la, ipiv);
  LAPACKE_dgetri(CblasRowMajor, la, M, la, ipiv);
  free(ipiv);

  // x(k+1) = x(k) + M^(-1)(b - Ax(k))
  while (res_norm > tol && counter < maxit) {

    cblas_dcopy(la, RHS, 1, work, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, -1.0, AB, lab, X, 1, 1.0, work, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, la, la, 1.0, M, la, work, 1, 1.0, X, 1);
    // Commented out for performance measurement
    // res_norm = cblas_dnrm2(la, work, 1);
    // res_norm /= cblas_dnrm2(la, X, 1);
    //printf("%zu %.16lf\n", counter, res_norm);
    counter++;
  }
  free(M);
  free(work);
}

void optimGaussSeidel(double *AB, double *RHS, double *X, int lab, int la, int ku, int kl,
                 double tol, int maxit) {
  double *work = malloc(sizeof(double) * la);

  size_t counter = 0;

  // This is suboptimal since M is a lower triangular
  // But I am too lazy to do this properly
  // The inverse of a tridiagonal matrix isn't tridiagonal
  // So we'll use a dense matrix to store the result (a lower triangular)
  double *M = malloc(sizeof(double) * la * la);

  // Build the lower diagonal matrix
  for (size_t i = 0; i < la; i++) {
    M[i * la + i] = AB[lab * i + 1];
  }

  for (size_t i = 0; i < la - 1; i++) {
    M[(i +1) * la + i] = AB[lab * i + 2];
  }

  // Inverse it
  int* ipiv = malloc(sizeof(int) * la);
  LAPACKE_dgetrf(CblasRowMajor, la, la, M, la, ipiv);
  LAPACKE_dgetri(CblasRowMajor, la, M, la, ipiv);
  free(ipiv);

  // x(k+1) = (x(k) - M^(-1) * Ax(k)) M^(-1)b
  while (counter < maxit) {

    cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1.0, AB, lab, X, 1, 0, work, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, la, la, -1.0, M, la, work, 1, 1.0, X, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, la, la, 1.0, M, la, RHS, 1, 1.0, X, 1);
    counter++;
  }
  free(M);
  free(work);
}


void richardsonWithAlpha(double *AB, double *RHS, double *X, double alpha_rich, int lab, int la, int ku, int kl,
                         double tol, int maxit) {
  double *work = malloc(sizeof(double) * la);

  double res_norm = 1;
  size_t counter = 0;

  // x(k+1) = x(k) + alpha(b - Ax(k))
  while (res_norm > tol && counter < maxit) {

    cblas_dcopy(la, RHS, 1, work, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, -1.0, AB, lab, X, 1, 1.0, work, 1);
    cblas_daxpy(la, alpha_rich, work, 1, X, 1);
    // Commented out for performance measurement
    //res_norm = cblas_dnrm2(la, work, 1);
    //res_norm /= cblas_dnrm2(la, X, 1);
    //printf("%zu %.16lf\n", counter, res_norm);
    counter++;
  }
  free(work);
}


