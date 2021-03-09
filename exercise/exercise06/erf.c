#include<gsl/gsl_integration.h>
#include<math.h>
#include<stdio.h>
#include<gsl/gsl_math.h>

double my_f(double x, void* params){
	double f = 2*exp(-pow(x,2))/sqrt(M_PI);
	return f;
}
	
double my_erf(double z){
	gsl_function F;
	F.function=&my_f;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc (1000);
	double a=0, epsabs=0, epsrel = 1e-6, result,error;
	size_t limit = 1000;
        gsl_integration_qags(&F,a,z,epsabs,epsrel,limit,w,&result,&error);
	return result;
}
int main(){
	double xmin = -5, xmax = 5;
	for(double x=xmin;x<=xmax;x+=0.1)
		printf("%10g %10g\n",x,my_erf(x));
return 0;
}
