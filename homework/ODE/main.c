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
	double b1_star = 1;
	double b2 = 1;
	for(int i = 0;i<n;i++){
		gsl_vector_set(yerr,i,h*(b2*gsl_vector_get(k2,i)-b1_star*gsl_vector_get(k0,i)));
		}
	gsl_vector_free(k0);
	gsl_vector_free(k1);
	gsl_vector_free(k2);
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

void fun(double t, gsl_vector* y, gsl_vector* dydt){
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
}

void SIR( double t, gsl_vector* y, gsl_vector* dydt){
	// SIR key-values
	// Note this model makes no assumption about quarantine as far as I know
	// The time between contact. I'm picking 7 day, because I have to pick something
	double Tc = 7;
	// On WebMD they say recovery is between 2 -14, so I'll pick 14
	double Tr = 14;
	double S = gsl_vector_get(y,0);
	double I = gsl_vector_get(y,1);
	double R = gsl_vector_get(y,2);
	double N = 5823000 -R;// -(I+ R); // THe number of resistant people is irrelavent
//	double R = gsl_Vector_get(y,2);
	gsl_vector_set(dydt,0,-I*S/(N*Tc));  // S
	gsl_vector_set(dydt,1,I*S/(N*Tc)-I/Tr);  // I
	gsl_vector_set(dydt,2,I/Tr);  // R
}
int main(){
	int n = 2;
	int m = 3;
	gsl_vector* y0   = gsl_vector_alloc(n);
	gsl_vector* yend = gsl_vector_alloc(n);
	gsl_vector* dydt = gsl_vector_alloc(n);
	gsl_vector* yerr = gsl_vector_alloc(n);
	gsl_vector* y01 = gsl_vector_alloc(m);
	gsl_vector* yend1 = gsl_vector_alloc(m);
	gsl_vector* dydt1 = gsl_vector_alloc(m);
	gsl_vector* yerr1 = gsl_vector_alloc(m);
	double dt = 0.1;
	double t0 = 0;
	double tend =4;
	gsl_vector_set(y0,0,0);
	gsl_vector_set(y0,1,1);
	double acc = 1e-2;
	double eps = 1e-2;
	fprintf(stderr,"\n\n");
	FILE* RESULTS = fopen("result.txt","w");
	fprintf(RESULTS,"\n # index 0 : The test of u''=-u\n");
	fclose(RESULTS); 
	driver(fun, t0, tend,dt,y0,yerr,acc,eps);
	fprintf(stderr,"\n\n");
	// The SIR initial conditions
	double N = 5823000; //Accounting to wikipedia the current population is around 5823000. 
	double I = 1; // we'll put the initial infected at 1
	gsl_vector_set(y01,0,N-I);
	gsl_vector_set(y01,1,I);
	gsl_vector_set(y01,2,0);
	FILE* RESULTS1 = fopen("result.txt","a");
	fprintf(RESULTS1,"\n\n# index 1 : SIR-model\n");
	fclose(RESULTS1); 
	double ts0 = 0;
	double tsend =365; 
	double dts = 1;
	double acc1 = 1e-4;
	double eps1 = 1e-4;
	driver(SIR, ts0, tsend,dts,y01,yerr1,acc1,eps1);
	gsl_vector_free(y0);
	gsl_vector_free(yend);
	gsl_vector_free(dydt);
	gsl_vector_free(yerr);
	gsl_vector_free(y01);
	gsl_vector_free(yend1);
	gsl_vector_free(dydt1);
	gsl_vector_free(yerr1);
	return 0;
}

