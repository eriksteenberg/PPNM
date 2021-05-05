#include<stdio.h>
#include"functions.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<gsl/gsl_blas.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#define RND (double)rand()/RAND_MAX


//double jacobi


void f(gsl_vector* xlist, gsl_vector* flist){ // this is f(x)= (x-1)**2
	int k=xlist->size;
	for(int i =0;i<k;i++){
		double x0 = gsl_vector_get(xlist,i);
		gsl_vector_set(flist,i,(x0-1)*(x0-1));
	}
}

void Rosen(gsl_vector* xlist, gsl_vector* flist){
	// We are to find the extremum, which we do by rooting the gradient, flist will be a list of gradients
	double x = gsl_vector_get(xlist,0); 
	double y = gsl_vector_get(xlist,1); 
	// f(x,y) = (1-x)**2 + 100*(y-x**2)**2
	//df/dx = 2*(1-x)*(-1) +2*100*(y-x**2)*(-2*x)
	//df/dy = 2*100*(y-x**2)*(1)
	//First dfdx
	gsl_vector_set(flist,0,2*(1-x)*(-1) + 2*100*(y-x*x)*(-2*x));
	gsl_vector_set(flist,1,2*100*(y-x*x));
}


void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){ 
	//for functions that takes a vector and returns a vector
	int n= x->size;
	gsl_matrix* jac = gsl_matrix_alloc(n,n);
	gsl_vector* dx = gsl_vector_alloc(n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* lambdax = gsl_vector_alloc(n);
	gsl_vector* lfx = gsl_vector_alloc(n);
	f(x,fx);
	int N = 0;
	while (N<1000){
		f(x,fx);
		jacobi(x,f,jac);
		gsl_vector_scale(fx,-1);
		GS_calculate(jac,fx,dx);
		N++;
		double lambda = 2;
		while(1){	
			lambda *=0.5;
			gsl_vector_scale(dx,lambda);
			f(x,fx);
			double norm_fx = gsl_blas_dnrm2(fx);
			gsl_vector_add(x,dx);
			f(x,lfx);//f(x+l*dx)
			double norm_lfx = gsl_blas_dnrm2(lfx);
			if(norm_lfx<(1-0.5*lambda)*norm_fx){//f(x+dx)>c*f(x)
			       //fprintf(stderr,"\n f(x+dx)<const.* f(x)\n");
			       break;}
			if(lambda>0.015625){
		//		fprintf(stderr,"\n Break because of lambda\n");
				break;}	
		}
	}
	gsl_matrix_free(jac);
	gsl_vector_free(dx);
	gsl_vector_free(fx);
	gsl_vector_free(lambdax);
	gsl_vector_free(lfx);
}

int main(){
	int n = 4;
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_matrix* jac = gsl_matrix_alloc(n,n);
	gsl_vector_set(x,0,1);	
	gsl_vector_set(x,1,-2);	
	gsl_vector_set(x,2,27);	
	gsl_vector_set(x,3,33);	
	printf("\n First we test our root-finder og f(x)=(x-1)**2)\nOur x vector\n");	
	show_vector(x);
	f(x,fx);
	printf("\nOur f(x) vector\n");
	show_vector(fx);
	/*jacobi(x,f,jac);
	show_matrix(jac);	
	GS_calculate(jac,fx,x);
	printf("vector x\n");
	show_vector(x);                        
	gsl_blas_dgemv(CblasNoTrans,1,jac,x,0,y);
	show_vector(y);
	void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){ //for funktions that take a vector and returns a vector
	*/
	double eps = 0.0001;
	newton(f,x,eps);
	printf("\n The root of the f(x) vector. It is a good sign that all the entries are 1\n");
	show_vector(x);
	printf("\n Now for the Test on the  Rosenbrock's valley function, f(x,y) = (1-x)**2 + 100*(y-x**2)**2,\n we search for the extremum points by applying our root finder on the derivatives of the function\n");
	int k = 2;
	gsl_vector* xy = gsl_vector_alloc(k);
	gsl_vector* fxy = gsl_vector_alloc(k);
	gsl_matrix* jac_1 = gsl_matrix_alloc(k,k);
	for(int i = 0;i<k;i++){
		gsl_vector_set(xy,i,RND);
	}
	printf("THese are randomly generated x0 and y0 values\n");
	show_vector(xy);
	newton(Rosen,xy,eps);
	printf("The return x and y values are\n");
	show_vector(xy);
	printf("Plucking them into the function, returns:\n");
	Rosen(xy,fxy);
	printf("\n When the derivates of a function are zero it is an extremum. So we have found what we wanted\n");
	show_vector(fxy);
	gsl_vector_free(fx);
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(xy);
	gsl_vector_free(fxy);
	gsl_matrix_free(jac);
	gsl_matrix_free(jac_1);
return 0;
}
