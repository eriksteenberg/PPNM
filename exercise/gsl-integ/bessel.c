#include<gsl/gsl_integration.h>
#include<math.h>
#include<stdio.h>
#include<gsl/gsl_math.h>

double f(double x, void* params){
	double n = *(double *) params;
	double z = *(double *) params;
	double f = cos(n*x-z*sin(x))/(M_PI);
        return f;
}

double mybessel(double z,double n){
	gsl_function F;
	F.function=&f;
	F.params=(void*)&n;
	F.params=(void*)&z;
	int limit=999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0,b=M_PI, acc=0, eps=1e-6,result,error;
	gsl_integration_qags(&F,a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}


int main(){
	double x = mybessel(5,0);
	double x1 = mybessel(10,0);
	printf("%g %g\n",x,x1);

return 0;
}
