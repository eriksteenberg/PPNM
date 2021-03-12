#include<gsl/gsl_linalg.h>
#include<stdio.h>
#include<math.h>
#include <gsl/gsl_integration.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#define RND (double)rand()/RAND_MAX
/*doublel interp(int n ,double* x ,double* y ,double z){ //lige nu er det vektor men vi kan også lave rækker med  double x[]
	assert(n>1 && z>=x[0] && z<=x[n−1]);
	int i =0,j=n−1;
	while( j−i>1){int m=(i+j)/2; if(z>x[m]) i=m; else j=m;}
	assert(x[i+1]>x[i]);
	return y[i]+(y[i+1]−y[i])/(x[i+1]-x[i])*(z-x[i])
}
*/
//struct params {gsl_vector x;gsl_vector y;};

double interp(double z,gsl_vector* x, gsl_vector* y){
	assert(z>gsl_vector_get(x,0) && z<=x->size);
	int i=0, j=x->size-1;
	while(j-1>1){int m=(i+j)/2; if(z>gsl_vector_get(x,m)) i=m; else j=m;}
	assert(gsl_vector_get(x,i+1)>gsl_vector_get(x,i));
	return gsl_vector_get(y,i)+(gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/(gsl_vector_get(x,i+1)-gsl_vector_get(x,i))*(z-gsl_vector_get(x,i));
}

double interp_integ(double z, gsl_vector* x, gsl_vector* y){
	double sum = 0;
	double dx = 1.0/10;
	for(double i=0.1;i<z;i+=dx){
		double dsum = dx * interp(i,x,y);
		sum+=dsum;
	}
	return sum;
}
/*void f(double z,void* arg){
	x=p->x;
	y=p->y;
        assert(z>gsl_vector_get(x,0) && z<=x->size);
        int i=0, j=x->size-1;
        while(j-1>1){int m=(i+j)/2; if(z>gsl_vector_get(x,m)) i=m; else j=m;}
        assert(gsl_vector_get(x,i+1)>gsl_vector_get(x,i));
	return gsl_vector_get(y,i)+(gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/(gsl_vector_get(x,i+1)-gsl_vector_get(x,i))*(z-gsl_vector_get(x,i));
}
double linterp_integ(double z,gsl_vector* x, gsl_vector* y){
	gsl_function F;
	F.function=&f;	
	int limit = 999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
        double a=0, b=1, acc=1e-6, eps=1e-6, result, error;
        gsl_integration_qags(&F, a, b, acc, eps, limit, w, &result, &error);
gsl_integration_workspace_free(w);
	return result;
}
*/
int main(){
	int n=2;
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	printf("Hej\n");
        for(int i=0; i< n; i++)
 		{
		double j = (double)i/10;
		gsl_vector_set(x,i,j);
		}
        for(int i=0; i< n; i++)
		{
		double bi=RND;
		gsl_vector_set(y,i,bi);
		}
	double sum1 = interp_integ(2.1,x,y);
	printf("%10g\n",sum1);
	double b = interp(1.4,x,y);
//	for(double i = 0.1;i<4.1;i+=1.0/8){
//		double sum = interp_integ(i,x,y);
//		printf("%g %g\n",i,sum);
//		}
gsl_vector_free(x);
gsl_vector_free(y);
return 0;
}
