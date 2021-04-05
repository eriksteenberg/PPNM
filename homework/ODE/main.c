#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#include"functions.c"
#define rnd (double)rand()/RAND_MAX

//void rkstep12(void f(int n, double x, double* yx, double* dydx), int n, double x, double* yx, double h, double* yh, double* dy)


void fun(double t, gsl_vector* y, gsl_vector* dydt){
	gsl_vector_set(dydt,0,1);
	gsl_vector_set(dydt,1,1);
//	gsl_vector_set(dydt,0,-gsl_vector_get(y,1));
//	gsl_vector_set(dydt,1,gsl_vector_get(y,0));
	}

void rkstep12(void f(double t, gsl_vector* y, gsl_vector* dydt), double t0, gsl_vector* y0, double h, gsl_vector* yend, gsl_vector* yerr){
	int n = y0->size;	//k0 is dydt
	gsl_vector* k0=gsl_vector_alloc(n);	//we generate dydt
	gsl_vector* k1=gsl_vector_alloc(n);	//we generate dydt
	f(t0,y0,k0);	//Now we calculate yend, by euler
	gsl_vector_memcpy(k1,y0);
	gsl_blas_daxpy(0.5*h,k0,k1);
	f(t0+0.5*h,k1,k1);
	for(int i=0;i<n;i++){
		double k1i = gsl_vector_get(k1,i);
		double y0i = gsl_vector_get(y0,i);
		gsl_vector_set(yend,i,y0i+h*k1i);
	}
	double b1_star = 1;
	double b2 = 1;
	for(int i = 0;i<n;i++){
		gsl_vector_set(yerr,i,h*(b2*gsl_vector_get(k1,i)-b1_star*gsl_vector_get(k0,i)));
		}
	gsl_vector_free(k0);
	gsl_vector_free(k1);
}

//void rkstep12(void f(double t, gsl_vector* y, gsl_vector* dydt), double t0, gsl_vector* y0, double h, gsl_vector* yend){
//	int n = y0->size;
//void fun1(double t, gsl_vector* y, gsl_vector* dydt){
//void fun2(double t, gsl_vector* y, gsl_vector* dydt){
int main(){
	int n = 2;
	gsl_vector* y0 = gsl_vector_alloc(n);
	gsl_vector* yend = gsl_vector_alloc(n);
	gsl_vector* dydt = gsl_vector_alloc(n);
	gsl_vector* yerr = gsl_vector_alloc(n);
	double dt = 0.1;
	double t0 = 0;
	for(int i=0;i<n;i++){
		gsl_vector_set(y0,i,rnd);
	}
	rkstep12(fun,t0,y0,dt,yend,yerr);
	printf("y0\n");
	show_vector(y0);
	printf("yend\n");
	show_vector(yend);
	printf("yerr\n");
	show_vector(yerr);
	gsl_vector_free(y0);
	gsl_vector_free(yend);
	gsl_vector_free(dydt);
	gsl_vector_free(yerr);
	return 0;
}
