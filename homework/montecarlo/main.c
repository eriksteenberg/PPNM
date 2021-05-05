#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<complex.h>
#define RND ((double)rand()/RAND_MAX)


void randomx(int dim, double *a, double *b, double *x){
	for (int i=0;i<dim;i++){
		x[i] = a[i] + RND*(b[i]-a[i]);
	}
}

// I made two versions of the monte carlo integration. plain this the one given in the exercise, the other is based upon table 1 in montecarlo.pdf. The one based upon the pdf was altered to be similar to the exercise on. THey ar eso similar that they are barely two sepereate functions, but I learned studd from implimenting both. SO theyere staying.

void plainmc(int dim,double *a, double *b, double f(int dim, double *x),int N, double* result, double* error){
	double V = 1;
	for(int i=0;i<dim;i++){
		V*=b[i]-a[i];
	}
	double sum=0, sum2=0, fx=0, x[dim];
	for(int i=0;i<N;i++){
		randomx(dim,a,b,x);	
		fx=f(dim,x);
		sum+=fx;
		sum2+=fx*fx;
	}
	double avr = sum/N, sigma2 = (sum2/N - avr*avr);
	//double sigma = sqrt(sum2/N - sum*sum/(N*N));	
	fprintf(stderr,"sigma2 is %g",sigma2);
	*result = (V) * avr;
	*error = sqrt(sigma2/N)*V;
}

complex plain(int dim, double f(int dim,double* x), double* a, double* b, int N){
	//The variables
	double V=1, sum=0, sum2=0, x[dim],mean=0,sigma=0;
	//The volume
	for(int i=0;i<dim;i++){
		V*=b[i]-a[i];
	}
	for(int i=0;i<N;i++){
		randomx(dim,a,b,x);
		double fx = f(dim,x);
		sum+=fx;
		sum2+=fx*fx;
	}
	mean = sum/N;
	sigma = sqrt(sum2/N-mean*mean);
	complex result = mean*V + I*sigma*V/sqrt(N);
	return result;
}
double f(int dim,double *x){
	double sum = 0;
	for(int i=0;i<dim;i++){
		sum+= x[i]*x[i];
	}
	//sum = sqrt(sum);
	return sum;
}

double f1(int dim,double *x){
	double sum = 1/(1 - cos(x[0]) * cos(x[1]) * cos(x[2]));
	return sum;
}

int main(){
	int k = 1;
	double a[k];
	double b[k];
	printf("The values of the lists that we are integrating over\n");
	for(int i=0;i<k;i++){
		a[i] = 0;
		b[i] = 5;
		printf("a[%i]=%g and b[%i]=%g\n",i,a[i],i,b[i]);
	}
	double error=0, result=0;
	int N = 1000000;
	plainmc(k,a,b,f,N,&result,&error);
	printf("mcplain integration of x**2 gives %g with an error of %g\n",result,error);
	//plain(int dim, double f(int dim,double* x), double* a, double* b, int N){
	
	// Now we calculate that crazy equation
	
	int n = 3;
	double a1[n];
	double b1[n];
	printf("The values of the lists that we are integrating over\n");
	for(int i=0;i<n;i++){
		a1[i] = 0;
		b1[i] = M_PI;
		printf("a[%i]=%g and b[%i]=%g\n",i,a1[i],i,b1[i]);
	}
	double error1=0, result1=0;
	int N1 = 1000000;
	plainmc(k,a1,b1,f1,N1,&result1,&error1);
	printf("montecarlo plain integration  of the crazy function gives %g with an error of %g\n",result1,error1);
	return 0;
}
