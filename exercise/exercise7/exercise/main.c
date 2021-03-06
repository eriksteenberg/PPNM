#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_sf_gamma.h>
#include<gsl/gsl_sf_erf.h>

double mygamma(double);

double myerf(double);

int main(){
	FILE* erf_out = fopen("erf.txt","w");
	FILE* gamma_out = fopen("gamma.txt","w");
	FILE* ln_gamma_out = fopen("ln.gamma.txt","w");

	double xmin = 0.1, xmax = 5;
	for(double x=xmin;x<xmax;x+=0.1){
		fprintf(ln_gamma_out,"%10g %10g %10g %10g \n",x,log(mygamma(x)),lgamma(x),log(gsl_sf_gamma(x)));
	}

//        double xmin =0.1, xmax = 5;
        for(double x=xmin;x<xmax;x+=0.1){
                fprintf(gamma_out,"%10g %10g %10g %10g \n", x, tgamma(x), gsl_sf_gamma(x), mygamma(x));
	}
  //      double xmin=-2, xmax=3.6;
        for(double x=xmin;x<xmax;x+=1.0/8){
                fprintf(erf_out,"%10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),myerf(x));}
fclose(erf_out);
fclose(gamma_out);
fclose(ln_gamma_out);
return 0;
}
