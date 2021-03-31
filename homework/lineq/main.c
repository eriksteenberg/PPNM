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

void show_matrix(gsl_matrix* A){
	int n = A->size1; //længde af vektorene
	int m = A->size2; // antal vektorer
//	printf("The matrix looks like this\n");
	printf("\n");
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			double x = gsl_matrix_get(A,i,j);
			printf("%9.3g ",x);
		}
		printf("\n");
	}
	printf("\n");
}

double cdot(gsl_vector* A, gsl_vector* B){
	assert(A->size==B->size);	
	double result=0;
	for(int i = 0;i<A->size;i++){
		double b = gsl_vector_get(A,i)*gsl_vector_get(B,i);
		result += b;	
	}
	return result;
}


void show_vector(gsl_vector* V){
	int n = V->size;
	printf("\n");
//	printf("The vector looks like\n");
	for(int i=0;i<n;i++){
		double x = gsl_vector_get(V,i);
		printf("%g\n",x);
	}
	printf("\n");
}
void GS_decomp(gsl_matrix* A,gsl_matrix* R){
	assert(A->size2 == R->size1);
	int N = A->size1;
	int M = A->size2;
//	gsl_matrix* Q = gsl_matrix_alloc(N,M);	

		gsl_vector* ai = gsl_vector_alloc(N);
		gsl_vector* aj = gsl_vector_alloc(N);
	for(int i=0;i<M;i++){
                gsl_matrix_get_col(ai,A,i);
		double Rii = gsl_blas_dnrm2(ai);
		gsl_matrix_set(R,i,i,Rii);
		gsl_vector_scale(ai,1/Rii);
		gsl_matrix_set_col(A,i,ai);
		for(int j=i+1;j<M;j++){
	  		gsl_matrix_get_col(aj,A,j);
			double Rij = cdot(ai,aj);
			gsl_blas_daxpy(-Rij,ai,aj);
			gsl_matrix_set(R,i,j,Rij);
			gsl_matrix_set_col(A,j,aj);
		}
	}
}
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	int n = x->size;
	gsl_vector* c = gsl_vector_alloc(n);
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,x);
/*	for(int i=0;i<n;i++){
		double Uii = gsl_matrix_get(R,i,i);
		gsl_vector_set(x,i,gsl_vector_get(c,i)/Uii);
		}
*/	for(int i=n-1;i>0;i--){
		double c1 = gsl_vector_get(x,i);
		for(int j=i+1;j<n;j++){
			c1 -= gsl_matrix_get(R,i,j)*gsl_vector_get(x,i);
		}
		gsl_vector_set(x,i,c1/gsl_matrix_get(R,i,i));
	}
	/*for(int i=n-1;i>0;i--){
		double c1 = gsl_vector_get(c,i);
		for(int j=i+1;j<n;j++){
			c1 -= gsl_matrix_get(R,i,j)*gsl_vector_get(x,i);
		}
		gsl_vector_set(x,i,c1/gsl_matrix_get(R,i,i));
	}*/
	gsl_vector_free(c);
}	

void matrix_multi(gsl_matrix* A, gsl_matrix* B,gsl_matrix* result){
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A,B,0.0,result);
}

void transposed_multi(gsl_matrix* A, gsl_matrix* result){
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0,A,A,0.0,result);
	

}
int main(){
	int N =6;
	int M =4;
	gsl_matrix* Q = gsl_matrix_alloc(N,M);	
	gsl_matrix* A = gsl_matrix_alloc(N,M);	
	gsl_matrix* B = gsl_matrix_alloc(M,M);	
	gsl_matrix* R = gsl_matrix_alloc(M,M);	
	gsl_matrix* ATA = gsl_matrix_alloc(N,M);	
	for(int i=0;i<N;i++){
    		for(int j=0;j<M;j++){
			gsl_matrix_set(Q,i,j,rnd);
		}
	}
//--------------------------------------//
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
//--------------------------------------------	
	int N1 = 2;
	printf("\n\n Now for GS_solve, we generate a vector b\n and a new square matrix A\n");
	//We didnt use V so we wil use it again.
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
	printf("\nAnd now We solve it\n");
	GS_solve(Q1,R1,V,x1);
	show_vector(x1);
	printf("checking Ax=b\n");
	gsl_vector* y1 = gsl_vector_alloc(N1);
	gsl_blas_dgemv(CblasNoTrans,1,T,x1,0,y1);
	printf("\n Ax is\n");
	show_vector(y1);
	printf("\n the vector b is \n");
	show_vector(V);
	printf("\nAs we can see, they fit fairly well \n");
	gsl_matrix_free(Q);
	gsl_matrix_free(Q1);
	gsl_matrix_free(R1);
	gsl_matrix_free(T);
	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_matrix_free(ATA);
	gsl_matrix_free(R);
	gsl_vector_free(V);
	gsl_vector_free(y1);
	gsl_vector_free(x1);
return 0;
}
