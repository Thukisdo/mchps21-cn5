#pragma once
#include <stdlib.h>
#include <stdio.h>

double* random_dvec(size_t m, double lower, double upper);
float* random_svec(size_t m, float lower, float upper);

void print_dvec(double* vec, size_t m);
void print_svec(float* vec, size_t m);