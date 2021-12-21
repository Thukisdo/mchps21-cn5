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

typedef struct {
  double *rhs;
  double *sol;
  double *exact_sol;
  double *x;
  int la;
} Poisson1DProblem;

Poisson1DProblem createPoisson1D(int la, double T0, double T1, int write_to_file);

void freePoisson1D(Poisson1DProblem p);

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

void gaussSeidel(double *AB, Poisson1DProblem p, int lab, int ku, int kl,
                 double tol, int maxit, FILE *out);

void richardsonWithAlpha(double *AB, Poisson1DProblem p, int lab, int ku, int kl,
                         double tol, int maxit, FILE *out);