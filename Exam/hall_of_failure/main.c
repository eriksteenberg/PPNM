#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX
#include"functions.h"
#include<time.h>
gsl_vector* d;
gsl_vector* u;
static int p;
void sum(gsl_vector* lambda,gsl_vector* eigen){
	int n = lambda->size;
	fprintf(stderr,"lambda = %i, u = \n",n);
	fprintf(stderr,"Sup\n");
	assert(lambda->size==d->size);
	fprintf(stderr,"Sup2\n");
	assert(lambda->size==u->size);
	fprintf(stderr,"Sup1\n");
	double uk2,uk,sum,dp,l;
	dp = gsl_vector_get(d,p);
	for(int i=0;i<n;i++){
		sum =0;
		l = gsl_vector_get(lambda,i);
		sum+=-(dp+2*gsl_vector_get(u,p)-l);
		for(int k=0;k<n;k++){
			uk = (gsl_vector_get(u,k));
			if(uk==0){
				uk2=0; //I had trouble with numbers being devided by zero which was fixed by this.
			}
			else{
				uk2 = uk * uk / (gsl_vector_get(d,k)-l);
			}
			if(k==p){uk2=0;} // 
			sum+=uk2;}
		gsl_vector_set(eigen,i,sum);}
}
void eigen(gsl_matrix* A, gsl_vector* x,int p){
	int n = A->size1;
	gsl_vector* e = gsl_vector_calloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(x,i,rnd);//10*rnd);
	}
	gsl_vector_set(e,p,1);
	printf("The initial lambda value guess\n");
	show_vector(x);
	gsl_blas_dger(1,e,u,A);
	gsl_blas_dger(1,u,e,A);
	printf("The e(p)-vector\n");
	show_vector(e);
	printf("The update Matrix\n");
	show_matrix(A);
	newton(sum,x,1e-8);
//void eigen(gsl_matrix* A, gsl_vector* x, gsl_vector* d, gsl_vector* u, int p)void sum(gsl_vector* lambda,gsl_vector* eigen,gsl_vector* d, gsl_vector* u, int p){
	
	gsl_vector_free(e);
}


int main(){
	int n = 5;
	p = 2;
	printf("We begin by first testing my implimentation of the Symmetric row/column update for a %ix%i-matrix for a random diagonal matrix annd update vector and e(p)=4(%i)\n",n,n,p);
	//Matrix being  updated
	gsl_matrix* A = gsl_matrix_calloc(n,n);
	// u-vector
	u = gsl_vector_calloc(n);
	// eigen valeu vector
	gsl_vector* x = gsl_vector_calloc(n);
	// a vector containing all the d_p values of the diagonal matrix
	d = gsl_vector_calloc(n);
	//a vector holding the lambda values
	gsl_vector* lambda = gsl_vector_calloc(n);
	//e(p) is introduced and freed in the eigen function
	double dpi =0;
	fprintf(stderr,"Hej\n");	
	for(int i =0;i<n;i++){
		dpi =rnd;
		gsl_matrix_set(A,i,i,dpi);
		gsl_vector_set(d,i,dpi);
		gsl_vector_set(lambda,i,rnd);
		gsl_vector_set(u,i,rnd);
		}
	gsl_vector_set(u,p,0);
	printf("The diagonal matrix is\n");
	show_matrix(A);
	printf("e(p)=e(%i) and the update vector is \n",p);
	show_vector(u);
	printf("Now we test my row/coloum update implimentation, or the eigen value calcualtor as I've called it, on the randomly generated matrix.\nThe function takes a diagonal matrix, an update vector, a vector to store the eigenvalues and a int p-value and then returns solutions to the eigenvalue equation (4.30) in the storing vector.\n");
	fprintf(stderr,"Hej\n");	
	eigen(A,x,p);
	fprintf(stderr,"Hej der\n");	

	printf("\nThe list of of values that eigen value calculator returned\n");
	show_vector(x);
	printf("We then pluck these values into eq. 4.30, just to check\n");
	sum(x,lambda);
	printf("The vector of the sums with the calculated eigenvalues\n");
	show_vector(lambda);
	printf("\n\n\nTo estimate the number of operations, the following was done:\nA couple of random matrices of fairly different orders were calculated and the time need to do the calculation was recorded. Then the time versus dimension of the matrix was plotted to compare.\n");
	gsl_vector_free(d);
	gsl_vector_free(u);
	gsl_matrix_free(A);
	gsl_vector_free(lambda);
	gsl_vector_free(x);
	return 0;
}
