#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#define RND ((double)rand()/RAND_MAX)


void randomx(int dim, double *a, double *b, double *x){
	for (int i=0;i<dim;i++){
		x[i] = a[i] + RND*(b[i]-a[i]);
	}
}
void plainmc(int dim,double *a, double *b, double f(double *x),int N, double* result, double* error){
	double V = 1;
	for(int i=0;i<dim;i++){
		V*=b[i]-a[i];
	}
	//fprintf(stderr,"V is %g\n",V);
	double sum=0, sum2=0, fx, x[dim];
	for(int i=0;i<N;i++){
		randomx(dim,a,b,x);
		fx=f(x);
		sum+=fx;
		sum2+=fx*fx;
	}
	double avr =sum/N, var = sum2/N - avr*avr;	
	//fprintf(stderr,"avr=%g var=%g sum2=%g sum=%g\n",avr,var,sum2,sum);
	//fprintf(stderr,"%g\n",avr*V);
	*result = avr*V;
	//fprintf(stderr,"result is %p\n",result);
	*error=sqrt(var/N)*V;
}

double f(double *x){
	return *x;
}
int main(){
	int k = 4;
	double a[k];
	double b[k];
	double x[k];
	for(int i=0;i<k;i++){
		a[i] = RND;
		b[i] = RND;
		printf("a[%i]=%g and b[%i]=%g\n",i,a[i],i,b[i]);
	}
	randomx(k,a,b,x);
	for(int i=0;i<k;i++){
	printf("the random x is %g\n",x[i]);
	}
	double error=0, result=0;
	int N = 100;
	plainmc(k,a,b,f,N,&result,&error);
	printf("mcplain gives %g with an error of %g\n",result,error);
return 0;
}
