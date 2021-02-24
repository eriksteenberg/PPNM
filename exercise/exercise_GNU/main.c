#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>

double erf(double);
int main(){
	double xmin=-3, xmax=3;
	for(double x=xmin;x<xmax;x+=1.0/8){
		printf("%10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),erf(x));
	}
return 0;	
}
