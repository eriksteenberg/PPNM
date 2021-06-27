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
//	int j = u->size;
//	fprintf(stderr,"lambda = %i, u = \n",n);
//	fprintf(stderr,"Sup\n");
	assert(lambda->size==d->size);
//	fprintf(stderr,"Sup2\n");
//	assert(lambda->size==u->size);
//	fprintf(stderr,"Sup1\n");
	double uk2,uk,sum,dp,l;
	dp = gsl_vector_get(d,p);
	for(int i=0;i<n;i++){
		sum =0;
		l = gsl_vector_get(lambda,i);
		sum+=-(dp+2*gsl_vector_get(u,p)-l);
	//	fprintf(stderr,"-(d_p - lambda) = %g",sum);
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
	gsl_blas_dger(1,e,u,A);
	gsl_blas_dger(1,u,e,A);
	newton(sum,x,1e-8);
	
	gsl_vector_free(e);
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
	double dpi =0;
	for(int i =0;i<n;i++){
		dpi =rnd;
		gsl_matrix_set(A,i,i,dpi);
		gsl_vector_set(d,i,dpi);
		gsl_vector_set(u,i,rnd);
		}
	gsl_vector_set(u,p,0);
	gettimeofday(&start, NULL);
	eigen(A,x,p);
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)*1e-6;
	printf("%i %.5f\n",n,runtime);
	gsl_vector_free(d);
   	gsl_vector_free(u);
	gsl_matrix_free(A);
	gsl_vector_free(x);
	}

int main(int arg, char** argv){
        n = (int)atof(argv[1]);
        p = 2;
	timeFile();
/*	clock_t start, end;
	double cpu;
	start = clock();
	eigen(A,x,p);
	end = clock();
       	
	cpu = ((double) (end - start)) / CLOCKS_PER_SEC;

	gsl_vector_free(d);
	gsl_vector_free(u);
	printf("%i %.25g\n",n,cpu);
	gsl_matrix_free(A);
	gsl_vector_free(lambda);
	gsl_vector_free(x);
	*/
	return 0;
}
