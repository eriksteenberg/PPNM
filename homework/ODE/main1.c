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
	//	gsl_vector_set(dydt,0,1);
	//	gsl_vector_set(dydt,1,1);
	gsl_vector_set(dydt,0,-gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,gsl_vector_get(y,0));
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
	//void rkstep12(void f(double t, gsl_vector* y, gsl_vector* dydt), double t0, gsl_vector* y0, double h, gsl_vector* yend, gsl_vector* yerr){

void driver(void f(double t, gsl_vector* y, gsl_vector* dydt), double b,double h, gsl_vector* xlist, gsl_matrix* listofy , double acc, double eps){
	int n = listofy->size1;
	gsl_vector* yerr = gsl_vector_alloc(n);
	gsl_vector* yb = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	int k = 0;
	double x,  err, normy,tol;
	//double s;
	double a = gsl_vector_get(xlist, 0);//tjek
	while(gsl_vector_get(xlist, k) < b){
		x = gsl_vector_get(xlist, k); 
		gsl_matrix_get_col( y, listofy, k);
		if(x + h > b) h = b - x;
       		rkstep12( f, x, y, h, yb, yerr);
	//	show(vector(yb);
//		s=0;
//		for(int i=0;i<n;i++)s+=(gsl_vector_get(yerr,i))*(gsl_vector_get(yerr,i));
		err = gsl_blas_dnrm2(yerr);//sqrt(s);
		fprintf(stderr,"err is %g\n",err);
//		s=0;
//		for(int i=0;i<n;i++)s+=(gsl_vector_get(yb,i))*(gsl_vector_get(yb,i));
		normy = gsl_blas_dnrm2(yb);//sqrt(s);
		fprintf(stderr,"norm of y is %g\n",normy);
		tol = (normy * eps + acc) * sqrt(h/(b-a));
		if(err<tol){
			k++;
			gsl_vector_set(xlist,k,x+h);
			gsl_matrix_set_col(listofy,k,yb);
		//	show_vector(xlist);
		}
		if(err>0)h*= pow(tol/err,0.25)*0.95;
		else h*=2;
	}
	gsl_vector_free(yerr);
	gsl_vector_free(yb);
	gsl_vector_free(y);
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
	double acc = 1e-6;
	double eps = 1e-6;
	double tend = 20;
	int m = (tend-t0)/dt;
	fprintf(stderr,"%i\n",m);
	gsl_vector* xlist = gsl_vector_alloc(m);
	gsl_matrix* ylist = gsl_matrix_alloc(n,m);
	gsl_matrix_set_col(ylist,0,y0);
	gsl_vector_set(xlist,0,t0);
	driver(fun, tend,dt,xlist, ylist ,acc,eps);
//	show_matrix(ylist);
	gsl_vector_free(y0);
	gsl_vector_free(yend);
	gsl_vector_free(dydt);
	gsl_vector_free(yerr);
	return 0;
}
