#include<gsl/gsl_linalg.h>
#include<stdio.h>
#include<math.h>
#include <gsl/gsl_integration.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#define RND (double)rand()/RAND_MAX

int binsearch(gsl_vector* x, double z){
	int i=0, j=x->size-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>gsl_vector_get(x,mid)) i=mid; else j=mid;
		}
	return i;
	}

double interp(gsl_vector* x, gsl_vector* y,double z){
	int x_len = x->size;
	assert(z>gsl_vector_get(x,0) && z<gsl_vector_get(x,x_len-1) && x_len >1);	
	int i = binsearch(x,z);
	double slope = (gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/(gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
	double value = gsl_vector_get(y,i) + slope*(z-gsl_vector_get(x,i));
	return value;
}

double linterp(gsl_vector* x, gsl_vector* y, double z){
}

double interp_integ(gsl_vector* x, gsl_vector* y,double z){
	int x_len = x->size;
	assert(z>gsl_vector_get(x,0) && z<gsl_vector_get(x,x_len-1) && x_len >1);
	int i = binsearch(x,z);			 
 	double slope = (gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/(gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
	double integ = gsl_vector_get(y,0)*(z - gsl_vector_get(x,0)) + 0.5*slope*pow((z - gsl_vector_get(x,0)),2);
	return integ;
}
double function(double z){
	double z0 = 5;
	double n = 1;
	//double resul = n*(2*z0*z - pow(z,2));
	double resul = n*(pow(z,2));
	return resul;
}

int main(){
	int n=10;
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
        for(int i=0; i< n; i++)
 		{
		double j = (double)i;
		gsl_vector_set(x,i,j);
		double bi=function(j);
		gsl_vector_set(y,i,bi);
		}
	double  n1=9;
	for(double i=0.1;i<n1;i+=0.1)
	{
		double j = (double)i;
		double b = interp(x,y,i);
		double integ = interp_integ(x,y,i);
		printf("%10g %10g %10g %10g\n ",i,b,function(i),integ); 
	}
	
	gsl_vector_free(x);
	gsl_vector_free(y);
return 0;
}
