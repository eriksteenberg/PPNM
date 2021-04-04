#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX

double norm(double x);

void show_matrix(gsl_matrix* A);

void show_vector(gsl_vector* V);

void Jtimes(gsl_matrix* A, int p, int q, double theta);

void timesJ(gsl_matrix* A, int p, int q, double theta);

int jacobi_diag(gsl_matrix* A, gsl_matrix* V);

