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
	double uk2,uk,sum,dp,l;
	dp = gsl_vector_get(d,p);
	for(int i=0;i<n;i++){
		sum =0;
		l = gsl_vector_get(lambda,i);
		sum+=-(dp+2*gsl_vector_get(u,p)-l);
		for(int k=0;k<n;k++){
			uk = (gsl_vector_get(u,k));
			uk2 = uk * uk / (gsl_vector_get(d,k)-l);
			if(k==p){uk2=0;} // 
			sum+=uk2;}
		gsl_vector_set(eigen,i,sum);}
}


void eigen(gsl_matrix* A, gsl_vector* x,int p,gsl_vector* e){
	int n = A->size1;
	for(int i=0;i<n;i++){
		gsl_vector_set(x,i,rnd);//10*rnd);
	}
	gsl_vector_set(e,p,1);
	gsl_blas_dger(1,e,u,A);
	gsl_blas_dger(1,u,e,A);
	newton(sum,x,1e-6);
//void eigen(gsl_matrix* A, gsl_vector* x, gsl_vector* d, gsl_vector* u, int p)void sum(gsl_vector* lambda,gsl_vector* eigen,gsl_vector* d, gsl_vector* u, int p){
	
}


int main(){
	int n = 5;
	p = 2;
	printf("We demontrate my implimentation of the Symmetric row/column update for a %ix%i-matrix \n",n,n);
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
	//e(p) 
	gsl_vector* e = gsl_vector_calloc(n);
	double dpi =0;
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
	printf("Now we test my row/coloum update implimentation, or the eigen value calcualtor as I've called it, on the randomly generated matrix.\nThe function takes a diagonal matrix, an update vector, a vector to store the eigenvalues and an int p-value and then returns solutions to the eigenvalue equation (4.30) in the storing vector.\n");
	eigen(A,x,p,e);
        printf("The update Matrix\n");
        show_matrix(A);
	printf("\nThe list of of values that eigen value calculator returned\n");
	show_vector(x);
	printf("We then pluck these values into eq. 4.30, just to check\n");
	sum(x,lambda);
	printf("The vector of the sums with the calculated eigenvalues\n");
	show_vector(lambda);
	printf("Comment here\nAs the reader can see some of the given values are the same, this is not neccessarily calls for concern. There are two reasons for this.\n#1: It is not a garuantee that the equation used to calculate the eigenvalues gives all the possible eigenvalue.\nAs a clear example of this, the equation does not work for diagonal equations, where the update vector is zero, as it only has one solution. THus it is reasonable to assume that not all eigenvalues are covered by this equation.\n#2:THe function that I used to solve the eigenvalue equation is my implimentation of the newton method for functions of vectors. Unfortunately this function appears to be biased towards values in the vacinity of the initial guess. THis is why the initial guess for the eigenvalues are randonly genereate, so as to avoid only getting the same eigenvalue in all the solutions.\n");
	
	gsl_vector_free(d);
	gsl_vector_free(u);
	gsl_matrix_free(A);
	gsl_vector_free(lambda);
	gsl_vector_free(x);
	return 0;
}
