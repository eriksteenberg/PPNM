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
					fprintf(stderr,"%9i ",0);
				}
				else{
					fprintf(stderr,"%9.3g ",x);
				}
		}
		printf("\n");
	}
	printf("\n");
}
void show_vector(gsl_vector* V){
	int n = V->size;
	printf("\n");
	for(int i=0;i<n;i++){
		fprintf(stderr,"%g\n",gsl_vector_get(V,i));}
	printf("\n");
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

void driver(void f(double t, gsl_vector* y, gsl_vector* dydt),double a,
	double b,double h, gsl_vector* y,gsl_vector* yerr, double acc,double eps){
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
			fprintf(RESULTS,"%3g %3g %3g\n",k,gsl_vector_get(y,0),gsl_vector_get(y,1));
		}
		if(err>0)h*= pow(tol/err,0.25)*0.95;
		else h*=2;
	}
	fprintf(RESULTS,"\n\n");
	fclose(RESULTS);
	gsl_vector_free(yb);
	}

void fun(double t, gsl_vector* y, gsl_vector* dydt){
//	gsl_vector_set(dydt,0,3*t*t);
//	double y1 = gsl_vector_get(y,1);	
//	double y0 = gsl_vector_get(y,0);	
//	gsl_vector_set(dydt,0,y1);//
//	gsl_vector_set(dydt,1,-y0);
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
}
int main(){
	int n = 2;
	gsl_vector* y0 = gsl_vector_alloc(n);
	gsl_vector* yend = gsl_vector_alloc(n);
	gsl_vector* dydt = gsl_vector_alloc(n);
	gsl_vector* yerr = gsl_vector_alloc(n);
	double dt = 0.1;
	double t0 = 0;
	gsl_vector_set(y0,0,0);
	gsl_vector_set(y0,1,1);
	double acc = 1e-1;
	double eps = 1e-1;
	double tend =5;
	FILE* RESULTS = fopen("result.txt","w");
	fprintf(RESULTS,"# index 0 : The test of x*2\n");
	fclose(RESULTS);
	driver(fun, t0, tend,dt,y0,yerr,acc,eps);
	gsl_vector_free(y0);
	gsl_vector_free(yend);
	gsl_vector_free(dydt);
	gsl_vector_free(yerr);
	return 0;
}

