/**********************************************/
/* lib_poisson1D.h                            */
/* Header for Numerical library developed to  */
/* solve 1D Poisson problem (Heat equation)   */
/**********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "blaslapack_headers.h"

void makeRowMajorGBand(double *AB, int lab, int la, int kv);

void makeColMajorGBand(double *AB, int lab, int la, int kv);

void makeColMajorGBIdentity(double *AB, int lab, int la, int kv);

void set_dense_RHS_DBC_1D(double *RHS, int *la, double *BC0, double *BC1);

void computeAnalyticalSolution(double *EX_SOL, double *X, int *la, double *BC0, double *BC1);

void setGridPoints(double *x, int *la);

void writeRowMajorGBand(double *AB, int *lab, int *la, char *filename);

void writeColMajorGBand(double *AB, int *lab, int *la, char *filename);

void writeVec(double *vec, int *la, char *filename);

void write_xy(double *vec, double *x, int *la, char *filename);

void computeEigenValues(double *eigval, int *la);

double computeMaxEigenValue(int la);

double computeMinEigenValue(int la);

double computeRichardsonOptAlpha(int la);

void gaussSeidel(double *AB, double *RHS, double *X, int lab, int la, int ku, int kl,
                 double tol, int maxit);

void optimGaussSeidel(double *AB, double *RHS, double *X, int lab, int la, int ku, int kl,
                 double tol, int maxit);

void richardsonWithAlpha(double *AB, double *RHS, double *X, double alpha_rich, int lab, int la, int ku, int kl,
                         double tol, int maxit);

void optimRichardsonWithAlpha(double *AB, double *RHS, double *X, double alpha_rich, int lab, int la, int ku, int kl,
                         double tol, int maxit);
