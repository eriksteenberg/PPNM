#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX
#define RND (double)rand()/RAND_MAX


double cdot(gsl_vector* A, gsl_vector* B);

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void GS_decomp(gsl_matrix* A,gsl_matrix* R);

void matrix_multi(gsl_matrix* A, gsl_matrix* B,gsl_matrix* result);

void transposed_multi(gsl_matrix* A, gsl_matrix* result);

void lsfit(int m, double f(int i, double x), gsl_vector* x, gsl_vector* y,gsl_vector* dy,gsl_vector* c, gsl_matrix* S);

void jacobi(gsl_vector* x, void f(gsl_vector* xlist, gsl_vector* flist),gsl_matrix* Jac);

void GS_calculate(gsl_matrix* A, gsl_vector* b, gsl_vector* x);

void jac_special(double f(gsl_vector* x),gsl_vector* x, gsl_matrix* jac);

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps);