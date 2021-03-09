#include<gsl/gsl_integration.h>
#include<math.h>
#include<stdio.h>
#include<gsl/gsl_math.h>

double f(double x, void* params){
	double f = log(x)/sqrt(x);
	return f;
	}
	
int main(){
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	double error, result;
	gsl_function F;
	F.function = &f;	
	int a = 0; //the starting point
	int b = 1;
	double epsabs =0, epsrel = 1e-7;
	size_t limit = 1000;
	gsl_integration_qags(&F,a,b,epsabs,epsrel,limit,w,&result,&error);
	printf("numerically integrating the function gives %10g\n",result);	
	return 0;
}
