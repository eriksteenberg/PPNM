#include<stdio.h>
#include"functions.h"
#include<gsl/gsl_vector.h>
//#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<gsl/gsl_blas.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#define RND (double)rand()/RAND_MAX


double f(gsl_vector* x){
	int m = x->size;
	double sum = 0;
	double xi = 0;
	double xi2 = 0;
	for(int i = 0;i<m;i++){
		xi = gsl_vector_get(x,i);
		xi2 = xi*xi;
		sum += xi2;
	}
	return sum;
}

double Rosen(gsl_vector* xi){
	double x = gsl_vector_get(xi,0);
	double y = gsl_vector_get(xi,1);
	return (1-x*x)*(1-x*x) + 100*(y-x*x)*(y-x*x);
}
double Himmel(gsl_vector* xi){
	double x = gsl_vector_get(xi,0);
	double y = gsl_vector_get(xi,1);
	return (x*x+y-11)*(x*x+y-11) + (x+y*y-7)*(x+y*y-7);	
}

double F(double E, double m, double gamma, double A){
	return A/((E-m)*(E-m)+gamma*gamma/4);
}
double SIGMA(double E,double m,double gamma, double A,double sigma,double dsigma){
	return (pow(F(E,m,gamma,A)-sigma,2)/(dsigma*dsigma));
}
double D(gsl_vector* list){
	double m = gsl_vector_get(list,0);
	double gamma = gsl_vector_get(list,1);
	double A = gsl_vector_get(list,2);
	double sum = 0;
	double E[] = {101,103,105,107,109,111,113,115,117,119,121,123,125,127,129,131,133,135,137,139,141,143,145,147,149,151,153,155,157,159};
	double sigma[] = {-0.25,-0.3,-0.15,-1.71,0.81,0.65,-0.91,0.91,-2.52,-1.01,2.01,4.83,4.58,1.26,1.01,-1.26,0.45,0.15,-0.91,-0.81,-1.41,1.36,0.5,-0.45,1.61,-2.21,-1.86,1.76,-0.5};
	double dsigma[] = {2,2,1.9,1.9,1.9,1.9,1.9,1.9,1.6,1.6,1.6,1.6,1.6,1.6,1.3,1.3,1.3,1.3,1.3,1.3,1.1,1.1,1.1,1.1,1.1,1.1,1.1,0.9,0.9,0.9};
	int M = 29;
	double Ei=0,sigmai=0,dsigmai=0;
	for(int i=0;i<M;i++){
		Ei = E[i];
		sigmai = sigma[i];
		dsigmai = dsigma[i];
		sum += SIGMA(Ei,m,gamma,A,sigmai,dsigmai);
	}
	return sum;
}


int main(){
	int n = 1;
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* grad = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(x,i,RND);
	}
	double eps = 1e-6;
	int steps = 1000;
	printf("Exercise A)\n We calculate the extremum of the Rosenbrock's valley and Himmelblau's functions\n");
	n = 2;
	gsl_vector* x0 = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(x0,i,RND);
	}
	printf("First we calculate Rosenbrock: \nf(x,y) =(1-x)**2 +100*(y-x**2)**2 \nTaking two random x y values ");
	show_vector(x0);
	printf("We then use the method\n");
	quasinewton(Rosen,x0,eps,steps);
	printf("The method tells us that the extremum is at ");
	show_vector(x0);
	printf("Rosenbrock at the extremum is f(%g,%g)=%g\n",gsl_vector_get(x0,0),gsl_vector_get(x0,1),Rosen(x0));
	for(int i=0;i<n;i++){
		gsl_vector_set(x0,i,RND);                                             }
	printf("Then we calculate Himmelblau:\n f(x,y)=(x**2+y-11)**2 +(x+y**2-7)**2\n Taking another two random x y values ");
	show_vector(x0);
	quasinewton(Himmel,x0,eps,steps);
	printf("The method tells us that the extremum is at");
	show_vector(x0);
	printf("Himmelblau at the extremum is f(%g,%g)=%g\n",gsl_vector_get(x0,0),gsl_vector_get(x0,1),Himmel(x0));
	printf("The gradient of the Himmel at x0 is");
	show_vector(grad);
	printf("And thus we have solved A)\n");
	printf("Exercise B)\n");
	n=3;
	double eps1 = 1e-3;
	steps = 10000;
	gsl_vector* H = gsl_vector_alloc(n);
	// 0 = m, 1 = gamma, 2 = A
	gsl_vector_set(H,0,125);
	gsl_vector_set(H,1,1);
	gsl_vector_set(H,2,10);
	quasinewton(D,H,eps1,steps);
	show_vector(H);
	double hej = D(H);
			
	printf("D(H)=%g\n",hej);
	gsl_vector_free(H);
	gsl_vector_free(x);
	gsl_vector_free(x0);
	gsl_vector_free(grad);
return 0;
}
