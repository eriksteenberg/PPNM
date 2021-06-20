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
	int N = 6;
	int M = 4;
	fprintf(stderr,"HEJ!!!!!!\n");
	gsl_matrix* Q = gsl_matrix_alloc(N,M);	
	gsl_matrix* A = gsl_matrix_alloc(N,M);	
	gsl_matrix* B = gsl_matrix_alloc(M,M);	
	gsl_matrix* R = gsl_matrix_alloc(M,M);	
	gsl_matrix* ATA = gsl_matrix_alloc(N,M);
//--------------------------------------//	
	for(int i=0;i<N;i++){
    		for(int j=0;j<M;j++){
			gsl_matrix_set(Q,i,j,rnd);
		}
	}
	gsl_matrix_memcpy(A,Q);
	printf("\n");
	printf("The matrix A\n");
	show_matrix(Q);	
	printf("\n");
	GS_decomp(Q,R);
	printf("The matrix Q\n");
	show_matrix(Q);	
	printf("\n");
	printf("As we can see R is upper triangular\n");
	show_matrix(R);	
	printf("\n\n\n");
	matrix_multi(Q,R,ATA);	
	printf("-------QR=A--------\n");
	show_matrix(ATA);
	printf("\n");
	transposed_multi(Q,B);
	printf("----QTQ----\n");
	show_matrix(B);	
	printf("\n");
	fprintf(stderr,"HEJ!!!!!!\n");
//--------------------------------------------	
	int N1 = 5;
	printf("\n\n Now for GS_solve, we generate a vector b\n and a new square matrix A\n");
	gsl_matrix* T = gsl_matrix_alloc(N1,N1);
	gsl_matrix* Q1 = gsl_matrix_alloc(N1,N1);
	gsl_matrix* R1 = gsl_matrix_alloc(N1,N1);
	gsl_vector* x1 = gsl_vector_alloc(N1);
	gsl_vector* V = gsl_vector_alloc(N1);
	for(int i=0;i<N1;i++){
		for(int j=0;j<N1;j++){
			gsl_matrix_set(T,i,j,rnd);
		}
		gsl_vector_set(V,i,rnd);	
	}
	gsl_matrix_memcpy(Q1,T);
	printf("\nTHe genereate A-matrix\n");
	show_matrix(Q1);
	printf("\nvector b is \n");
	show_vector(V);
	GS_decomp(Q1,R1);
	printf("\nThe Q-matrix\n");
	show_matrix(Q1);
	printf("\nThe R-matrix\n");
	show_matrix(R1);
	printf("\nAnd now We solve QRx=b for x, whcih gives\n");
	GS_solve(Q1,R1,V,x1);
	show_vector(x1);
	printf("checking Ax=b\n");
	gsl_vector* y1 = gsl_vector_alloc(N1);
	gsl_blas_dgemv(CblasNoTrans,1,T,x1,0,y1);
	printf("\n Ax is\n");
	show_vector(y1);
	printf("\n the vector b is \n");
	show_vector(V);
	printf("\nAs we can see, they fit\n");
	/*
	THis is here in case I want to check what Rx and Q^Tv. I think.
	gsl_vector* y2 = gsl_vector_alloc(N1);
	gsl_vector* y3 = gsl_vector_alloc(N1);
	gsl_blas_dgemv(CblasNoTrans,1,R1,x1,0,y2);
	gsl_blas_dgemv(CblasTrans,1,Q1,V,0,y3);
	
	show_vector(y2);
	show_vector(y3);

	*/
	gsl_matrix* P = gsl_matrix_alloc(N1,N1);
	gsl_matrix* I = gsl_matrix_alloc(N1,N1);
	gsl_matrix* Pinv = gsl_matrix_alloc(N1,N1);
	for(int i=0;i<N1;i++){
		for(int j=0;j<N1;j++){
			gsl_matrix_set(P,i,j,rnd);					
		}
	}
	printf("\nWe now test my GS_invert void function\nThe generated matrix\n");
	GS_inv(P, Pinv);
	show_matrix(P);
	printf("\nTHe inverted matrix\n");
	show_matrix(Pinv);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,P,Pinv,0,I);
	printf("\nThe matrix product of the two matrixes gives:\n");
	show_matrix(I);
	printf("\nWhich tells us that it works\n");
	gsl_matrix_free(Q);
	gsl_matrix_free(Q1);
	gsl_matrix_free(R1);
	gsl_matrix_free(T);
	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_matrix_free(P);
	gsl_matrix_free(I);
	gsl_matrix_free(Pinv);
	gsl_matrix_free(ATA);
	gsl_matrix_free(R);
	gsl_vector_free(V);
	gsl_vector_free(y1);
	//gsl_vector_free(y2);
	//gsl_vector_free(y3);
	gsl_vector_free(x1);
return 0;
}


