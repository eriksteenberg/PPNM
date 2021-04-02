#include"functions.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<gsl/gsl_blas.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX


int main(){
	int n=3;
	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* Anew = gsl_matrix_alloc(n,n);
	gsl_matrix* V = gsl_matrix_alloc(n,n);
	gsl_matrix* res = gsl_matrix_alloc(n,n);
	gsl_matrix* D = gsl_matrix_alloc(n,n);
	//gsl_matrix* V = gsl_matrix_alloc(n,n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			gsl_matrix_set(A,i,j,rnd);
		}
	}
	printf("before jacobi\n");
	show_matrix(A);
	gsl_matrix_memcpy(Anew,A);
	jacobi_diag(Anew,V);
/*	rintf("\nAfter jacobi\nA is\n");
	show_matrix(Anew);
	printf("\n V is\n");
	show_matrix(V);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,A,V,0,res);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,res,0,D);
	printf("\nVTAV\n");
	show_matrix(D);
*/	gsl_matrix_free(A);	
	gsl_matrix_free(Anew);	
	gsl_matrix_free(V);	
	gsl_matrix_free(res);	
	gsl_matrix_free(D);	
return 0;
}
