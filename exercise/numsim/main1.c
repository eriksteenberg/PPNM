#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double f (double x, void * params) {
	double f log(x) / sqrt(x);
	return f;
}

int main (void){
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double result, error;
	gsl_function F;
	F.function = &f;
	gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,w, &result, &error);
gsl_integration_workspace_free (w);
return 0;
}
