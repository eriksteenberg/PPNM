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
#include<sys/time.h>
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
	newton(sum,x,1e-3);
}
static int n;

void timeFile(void){
	//Creating object
	struct timeval start, end;
	double runtime;
	gsl_matrix* A = gsl_matrix_calloc(n,n);
	u = gsl_vector_calloc(n);
	gsl_vector* x = gsl_vector_calloc(n);
	d = gsl_vector_calloc(n);
	gsl_vector* e = gsl_vector_calloc(n);
	double dpi =0;
	for(int i =0;i<n;i++){
		dpi =rnd;
		gsl_matrix_set(A,i,i,dpi);
		gsl_vector_set(d,i,dpi);
		gsl_vector_set(u,i,rnd);
		}
	gsl_vector_set(u,p,0);
	gettimeofday(&start, NULL);
	eigen(A,x,p,e);
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)*1e-6;
	printf("%i %.5f %i\n",n,runtime,p);
	gsl_vector_free(d);
   	gsl_vector_free(u);
	gsl_matrix_free(A);
	gsl_vector_free(x);
	gsl_vector_free(e);
	}

int main(int arg, char** argv){
        n = (int)atof(argv[1]);
	if(n%2==0){p = (int)n*0.5;}
	else{p = (int)(1.+n)*0.5;}
	timeFile();
	return 0;
}
