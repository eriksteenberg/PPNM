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
	printf("Part A of this exercise, involves creating a \njacobian diagonalization operator and checking that VTV = 1,\nVTAV=D and VDVT =A\nAll of which is shown below\n");
	int n=5;
	gsl_matrix* A = gsl_matrix_calloc(n,n);
	gsl_matrix* V = gsl_matrix_calloc(n,n);
	gsl_matrix* res0 = gsl_matrix_calloc(n,n);
	gsl_matrix* res1 = gsl_matrix_calloc(n,n);
	gsl_matrix* res2 = gsl_matrix_calloc(n,n);
	gsl_matrix* D = gsl_matrix_calloc(n,n);
	gsl_matrix_set_identity(V); //seeting the Idnetity matrix
	//generating random matrix, that is equal to ints own transpose
	for(int i=0;i<n;i++){
	/*	for(int j=0;j<n;j++){
			if(i==j){
				gsl_matrix_set(A,i,j,rnd);
			}
			else{
				double x = rnd;
				gsl_matrix_set(A,i,j,x);
				gsl_matrix_set(A,j,i,x);
			}
		}
	*/
		gsl_matrix_set(A,i,i,rnd);
	}
	printf("before jacobi, A is\n");
	show_matrix(A);
	gsl_matrix_memcpy(D,A);
	//jacobi diagonalization
	jacobi_diag(D,V);
	printf("\nAfter jacobi D is\n");
	show_matrix(D);
	printf("\nAnd V is\n");
	show_matrix(V);

	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,V,0,res0);
	printf("VTV is\n");
	show_matrix(res0);
	
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,A,V,0,res0);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,res0,0,res1);
	printf("\nVTAV is\n");
	show_matrix(res1);
	printf("which is equal to D\n");
	show_matrix(D);
	
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,D,V,0,res0);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,res0,0,res2);
	printf("\nVTDV is\n");
	show_matrix(res2);
	printf("Compared with A, which is \n");
	show_matrix(A);


	gsl_matrix_free(res0);	
	gsl_matrix_free(res1);	
	gsl_matrix_free(res2);	
	gsl_matrix_free(A);	
	gsl_matrix_free(V);	
	gsl_matrix_free(D);	
return 0;
}
