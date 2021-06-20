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
	// --- setting up the matrix representation of H --- //
	int n = 50;
	double s = 1.0/(n+1);
	gsl_matrix* H = gsl_matrix_alloc(n,n);
	for(int i=0;i<n-1;i++){
		gsl_matrix_set(H,i,i,-2);
		gsl_matrix_set(H,i,i+1,1);
		gsl_matrix_set(H,i+1,i,1);
	}
	gsl_matrix_set(H,n-1,n-1,-2);
	gsl_matrix_scale(H,-1/(s*s));
	show_matrix(H);
	gsl_matrix* V = gsl_matrix_alloc(n,n);
	gsl_matrix_set_identity(V);
	jacobi_diag(H,V);
	printf("\nV\n");
	show_matrix(V);
	printf("\nDiagonalized H\n");
	show_matrix(H);
	for(int k=0;k<n/3;k++){
		double exact = M_PI*M_PI*(k+1)*(k+1);
		double calculated = gsl_matrix_get(H,k,k);
		printf("%i %g %g\n",k,calculated,exact);
	}
	FILE* eigen = fopen("eigen.txt","w");
//	fprint
	for(int k=0;k<3;k++){
		fprintf(eigen,"# index %i:%i order eigenfunc\n\n",k,k);
		fprintf(eigen,"0 0\n");
		for(int i=0;i<n;i++){
			fprintf(eigen,"%g %g\n",(i+1.0)/(n+1), gsl_matrix_get(V,i,k));
		}
		fprintf(eigen,"1 0\n\n\n");
	}
	fclose(eigen);
	gsl_matrix_free(H);
	gsl_matrix_free(V);
	
	return 0;
}
