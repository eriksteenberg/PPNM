#include<stdio.h>
#include<math.h>
//#include<gsl_integration.h>
#include<gsl/gsl_integration.h>
double f(double x, void* params){
        double func = (log(x)/sqrt(x));
return func;
}

double ing() {
	gsl_function F;
	F.function = &f;
	int limit = 999;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
        double a=0, b=1, acc=1e-6, eps=1e-6, result, error;
        gsl_integration_qags(&F, a, b, acc, eps, limit, w, &result, &error);
gsl_integration_workspace_free(w);
return result;
}
double f_erf(double x, void* params){
	double f = 2/sqrt(M_PI)*exp(-pow(x,2));
	return f;
}
double erf(double z){
	gsl_function F;
	F.function =&f_erf;
	int limit = 999;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
        double a=0, acc=1e-6, eps=1e-6, result, error;
        gsl_integration_qags(&F, a, z, acc, eps, limit, w, &result, &error);
gsl_integration_workspace_free(w);
return result;
}


int main(){
	printf("Exercise A: Integration the function gives %g and when I use an online calculator I get the same result\n",ing());
	FILE* erf_plot= fopen("erf.data.txt","w");
	for(double z=-5;z<=5;z+=0.1){
		fprintf(erf_plot,"%10g %10g\n",z,erf(z));
	}
fclose(erf_plot);
return 0;
}
