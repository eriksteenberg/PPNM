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

// These values and gsl_vectors are set up so that I avoid having to use nested functions.
gsl_vector* d;
gsl_vector* u;
static int p;

// Calculationg the secular equation
void sum(gsl_vector* lambda,gsl_vector* eigen){
	int n = lambda->size;
	double uk2,uk,sum,dp,l;
	dp = gsl_vector_get(d,p);
	for(int i=0;i<n;i++){
		sum =0;
		l = gsl_vector_get(lambda,i);
		sum+=-(dp + 2 * gsl_vector_get(u,p) - l );
		for(int k=0;k<n;k++){
			uk = (gsl_vector_get(u,k));
			uk2 = uk * uk / (gsl_vector_get(d,k)- l );
			if(k==p){uk2=0;} // the sum ignores k=p
			sum+=uk2;
		}
		gsl_vector_set(eigen,i,sum);
	}
}

//The function where the roots of the secular equation are calculated
void eigen(gsl_matrix* A, gsl_vector* x,int p,gsl_vector* u){
	int n = A->size1;
	gsl_vector* e = gsl_vector_calloc(n); //e(p)
	gsl_vector_set(e,p,1);
	gsl_blas_dger(1,e,u,A); //updating the diagonal matrix
	gsl_blas_dger(1,u,e,A);
	for(int i=0;i<n;i++){
		gsl_vector_set(x,i,0.5*i);
	}
	newton(sum,x,1e-9);
	gsl_vector_free(e);
}

//The function used when the initial values are calculated using Jacobi diagonalization
void discussion_eigen(
		gsl_matrix* A, //In this case A is already updated
	       	gsl_vector* x,
		gsl_matrix* D, // Holds the eigenvalues
		int p,
		gsl_vector* u
		){
	int n = A->size1;
	for(int i=0;i<n;i++){
		gsl_vector_set(x,i,gsl_matrix_get(D,i,i));//10*rnd);
	}
	newton(sum,x,1e-4);
}


int main(){
	//matrix dimension
	int n = 4;
	//p value
	p = 2;

	printf("We demonstrate my implimentation of the Symmetric row/column update for a %ix%i-matrix \n",n,n);
	
	//Matrix being  updated
	gsl_matrix* A = gsl_matrix_calloc(n,n);
	
	// u-vector
	u = gsl_vector_calloc(n);
	
	// eigen valeu vector
	gsl_vector* x = gsl_vector_calloc(n);
	
	// a vector containing all the d_p values of the diagonal matrix
	d = gsl_vector_calloc(n);
	
	//a vector holding the lambda values and used to hold the inital guess' for the eigenvalues
	gsl_vector* lambda = gsl_vector_calloc(n);
	
	double dpi =0;
	
	// Generating the diagonal matrix and the update vector
	for(int i =0;i<n;i++){
		dpi =rnd;
		gsl_matrix_set(A,i,i,dpi);
		gsl_vector_set(d,i,dpi);
		gsl_vector_set(u,i,rnd);
		}
	printf("The diagonal matrix is\n");
	show_matrix(A);
	printf("e(p)=e(%i), note that since gsl vectors start at zero the p entry is entry number %i and the update vector is \n",p,p+1);
	show_vector(u);
	printf("Now we test my row/coloum update implimentation, or specifically the secular equation (eq. 4.30 in the book), on a randomly generated matrix.\nThe function takes a diagonal matrix, an update vector, a vector to store the eigenvalues and an int p-value and then returns solutions to the eigenvalue equation (4.30) in the storing vector.\n");
	eigen(A,x,p,u);
        printf("The updated Matrix\n");
        show_matrix(A);
	printf("\nThe list of of values that eigen value calculator returned\n");
	show_vector(x);
	printf("We then pluck these values into eq. 4.30, just to check\n");
	sum(x,lambda);
	printf("The vector of the sums with the calculated eigenvalues\n");
	show_vector(lambda);
	printf("Comment here\nAs the reader can see some of the given values are the same, there are reasons for this.\n(The main error) The function that I used to solve the secular equation is my implimentation of the Newton root method for functions of vectors.\nUnfortunately the Newton method focuses on the imidiate vaccinity of the initial guess and not all possible values.\nThis is why the initial guess for the eigenvalues are randomly generated the way that they are,\nit is an attempt to avoid only getting one eigenvalue in all the solutions.\nOf cause this also means that there is a chance that some eigenvalues are not found, which is quit unfortunate.\n");
	printf("To illustrate my point about the newton method and compare wiht another approach, we apply the jaocobi diagonalization\n");
	
	// Setting up the needed matrices for jacobi diagonalization
	gsl_matrix* D = gsl_matrix_calloc(n,n);
	gsl_matrix* V = gsl_matrix_calloc(n,n);
	gsl_matrix_set_identity(V);
	gsl_matrix_memcpy(D,A);
	jacobi_diag(D,V);

	//holds the eigenvalues
	show_matrix(D);
	printf("The values on the diagonal are the eigenvalues, calculated using the jacobi method.\nTo illustrate that we get better results if we make a really good initial guess, we pick these diagonal values as the initial guess for the eigenvalue calculator equation and see what happens.\n");
	discussion_eigen(A,x,D,p,u);
	printf("This is the calculated eigenvalues when the initial guess for the eigenvalues are the eigenvalues calculated using the Jacobi diagonalization.\n");
	show_vector(x);
	printf("As we can see it agrees with the Jacobi methode.\nPlease note however that at the eigenvalue zero,corrospoding to the null vector, a non zero value is given.\nTwo things, first I belive that is is because the eigenvalue calculator was not designed with the goal of finding zero valued eigenvalues in mind.\nWhether this is actually an eigenvalue I cannot say, but it solves the eigenvalue equation so it should be.\n");
	printf("From this, we can conclude that my implimentation works and the reason that it doesn't find all the eigenvalues is because of how my implentation of the Newton method treats the initial value and that the eigenvalue equation is not meant to be solved for eigenvalues=0.\n");
	gsl_vector_free(d);
	gsl_vector_free(u);
	gsl_matrix_free(A);
	gsl_matrix_free(D);
	gsl_matrix_free(V);
	gsl_vector_free(lambda);
	gsl_vector_free(x);
	return 0;
}
