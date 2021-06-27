#include"functions.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX

double norm(double x){
	double result = sqrt(pow(x,2));
	return result;
}
void show_matrix(gsl_matrix* A){
	int n = A->size1; 
	int m = A->size2; 
	fprintf(stderr,"\n");
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			double x = gsl_matrix_get(A,i,j);
				if(norm(x)<10e-10){	// I added this part because seeing 10-16 5 times, hurts my eyes
					fprintf(stderr,"%9i ",0);}
				else{fprintf(stderr,"%9.3g ",x);}
		}
		printf("\n");
	}
	printf("\n");
}
void show_vector(gsl_vector* V){
	int n = V->size;printf("\n");
	for(int i=0;i<n;i++){
		fprintf(stderr,"%g\n",gsl_vector_get(V,i));}
	printf("\n");
}

void rkstep12(void f(double t, gsl_vector* y, gsl_vector* dydt), double t0, gsl_vector* y0, double h, gsl_vector* yend, gsl_vector* yerr){
	int n = y0->size;	//k0 is dydt
	gsl_vector* k0=gsl_vector_alloc(n);	//we generate dydt
	gsl_vector* k1=gsl_vector_alloc(n);	//we generate dydt
	gsl_vector* k2=gsl_vector_alloc(n);	//we generate dydt
	f(t0,y0,k0);	//Now we calculate yend, by euler
	gsl_vector_memcpy(k1,y0);
	gsl_blas_daxpy(0.5*h,k0,k1);
	f(t0+0.5*h,k1,k2);
	for(int i=0;i<n;i++){
		double k1i = gsl_vector_get(k2,i);
		double y0i = gsl_vector_get(y0,i);
		gsl_vector_set(yend,i,y0i+h*k1i);
	}
	double b1_star = 1;double b2 = 1;
	for(int i = 0;i<n;i++){
		gsl_vector_set(yerr,i,h*(b2*gsl_vector_get(k2,i)-b1_star*gsl_vector_get(k0,i)));
		}
	gsl_vector_free(k0);gsl_vector_free(k1);gsl_vector_free(k2);
}

void driver(void f(double t, gsl_vector* y, gsl_vector* dydt),double a, double b,double h, gsl_vector* y,gsl_vector* yerr, double acc,double eps){
	FILE* RESULTS = fopen("result.txt","a");
	int n = y->size;
	double err, normy,tol,k=a;
	gsl_vector* yb = gsl_vector_alloc(n);
	while(k<b){
		if(k+h>b)h=b-k;
       		rkstep12(f,k,y,h,yb,yerr);
		err = gsl_blas_dnrm2(yerr);
		normy = gsl_blas_dnrm2(yb);
		tol = (normy * eps + acc) * sqrt(h/(b-a));
		if(err<tol){
			k += h;	
			gsl_vector_memcpy(y,yb);
			fprintf(RESULTS,"%11g ",k);
			for(int i=0;i<n;i++){
				fprintf(RESULTS,"%11g ",(gsl_vector_get(y,i)));
			}
			fprintf(RESULTS,"\n");

		}
		if(err>0)h*= pow(tol/err,0.25)*0.95;
		else h*=2;
	}
	fprintf(RESULTS,"\n\n");
	fclose(RESULTS);
	gsl_vector_free(yb);
	}

// The differential equation that defines f
static double E;
void f(double r, gsl_vector* y, gsl_vector* dydt){
	gsl_vector_set(dydt,0,gsl_vector_get(y,0));
	gsl_vector_set(dydt,1, -2*(E - 1/r)*gsl_vector_get(y,0));
}

double Fe(double e, double r){
	double rmin = 1e-3;
	E=e;
	if(r<rmin)return r-r*r;
	gsl_vector* y = gsl_vector_alloc(2);
	gsl_vector_set(y,0,rmin-rmin*rmin);
	gsl_vector_set(y,1,1-2*rmin);
	gsl_vector* yerr = gsl_vector_alloc(2);
	driver(f,rmin,r,0.01,y,yerr,1e-1,1e-1);
	double y1 = gsl_vector_get(y,0);
	gsl_vector_free(y);
	gsl_vector_free(yerr);
	return y1;
}
static double rmax;

void M(gsl_vector* x, /* x*/
	gsl_vector* M_result /* fx */
	){
	double E1 = gsl_vector_get(x,0);
	assert(E1<0);
	gsl_vector_set(M_result,0,Fe(E1,rmax));
}


int main(){
	rmax = 8;
	double A = Fe(1,1);
	printf("%g",A);
	gsl_vector* M_result = gsl_vector_alloc(1);
	gsl_vector* x = gsl_vector_alloc(1);
	gsl_vector_set(x,0,-0.5);
	newton(M,x,1e-2);
	M(x,M_result);
	E = gsl_vector_get(x,0);
	fprintf(stderr,"e=%g",E);
//	show_vector(M_result);
	gsl_vector_free(M_result);
	gsl_vector_free(x);

	return 0;
}
