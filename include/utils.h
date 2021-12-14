#pragma once
#include <stdlib.h>
#include <stdio.h>

double* rand_dvec(size_t m, double lower, double upper);
double* rand_dmat(size_t m, size_t n, double lower, double upper);

void print_dvec(double* vec, size_t m);
void print_dmat(double* mat, size_t n, size_t m);