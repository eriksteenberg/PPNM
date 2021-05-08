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

double Higgs(double E, double m, double gamma, double A){
	return A/((E-m)*(E-m)+gamma*gamma/4);
}
double SIGMA(double E,double m,double gamma, double A,double sigma,double dsigma){
	return (pow(Higgs(E,m,gamma,A)-sigma,2)/(dsigma*dsigma));
}
gsl_vector* E;
gsl_vector* sigma;
gsl_vector* dsigma;

double D(gsl_vector* list){
	double m = gsl_vector_get(list,0);
	double gamma = gsl_vector_get(list,1);
	double A = gsl_vector_get(list,2);
	double sum = 0;
	double Ei=0,sigmai=0,dsigmai=0;
	for(int i=0;i<30;i++){
		Ei = gsl_vector_get(E,i);
		sigmai = gsl_vector_get(sigma,i);
		dsigmai = gsl_vector_get(dsigma,i);
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
	double eps = 1e-3;
	int steps = 1000;
	n = 2;
	gsl_vector* x0 = gsl_vector_alloc(n);
	printf("Exercise A)\n We calculate the extremum of the Rosenbrock's valley and Himmelblau's functions\n");
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
	/////////////////////////////////////////////
	printf("Exercise B)\n");
	n=30;
	E = gsl_vector_alloc(n);
	sigma = gsl_vector_alloc(n);
	dsigma = gsl_vector_alloc(n);
	FILE * dataset;
	dataset = fopen("data.txt","r");
	double Ei=0,sigmai=0,dsigmai=0;
	//fprintf(stderr,"importer data\n");
	for(int i=0;i<30;i++){
		int dummy = fscanf(dataset,"%lg %lg %lg",&Ei,&sigmai,&dsigmai);
		//fprintf(stderr,"%g %g %g\n",Ei,sigmai,dsigmai);
		gsl_vector_set(sigma,i,sigmai);
		gsl_vector_set(E,i,Ei);
		gsl_vector_set(dsigma,i,dsigmai);
	}
	fclose(dataset);
	double eps1 = 1e-6;
	steps = 10000;
	n = 3;
	gsl_vector* H = gsl_vector_alloc(n);
	// 0 = m, 1 = gamma, 2 = A
	gsl_vector_set(H,0,125);
	gsl_vector_set(H,1,5);
	gsl_vector_set(H,2,10);
	quasinewton(D,H,eps1,steps);
	printf("Our quasinewton parameters\n");
	show_vector(H);
	printf("mass is %.9g\n",gsl_vector_get(H,0));
	//double Higgs(double E, double m, double gamma, double A){
	double m = gsl_vector_get(H,0);
	double GAMMA = gsl_vector_get(H,1);
	double A = gsl_vector_get(H,2);
	double hej = Higgs(gsl_vector_get(E,1),m,GAMMA,A);
	printf("Higgs(H)=%g\n",hej);
	printf("D=%g\n",D(H));
	//show_vector(grad);
	FILE * DATA = fopen("plot.txt","w");
	fprintf(DATA,"# - data\n");
	for (double i =101;i<159;i+=0.1){
		//double 
		fprintf(DATA,"%g %g\n",i,Higgs(i,m,GAMMA,A));	
	}
	fclose(DATA);
	gsl_vector_free(H);
	gsl_vector_free(x);
	gsl_vector_free(x0);
	gsl_vector_free(grad);
return 0;
}
