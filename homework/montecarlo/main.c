#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<complex.h>
#include<assert.h>
#define RND ((double)rand()/RAND_MAX)


void randomx(int dim, double *a, double *b, double *x){
	for (int i=0;i<dim;i++){
		x[i] = a[i] + RND*(b[i]-a[i]);
	}
}

// I made two versions of the monte carlo integration. plain this the one given in the exercise, the other is based upon table 1 in montecarlo.pdf. The one based upon the pdf was altered to be similar to the exercise on. THey are so similar that they are barely two sepereate functions, but I learned stuff from implimenting both, so theye are staying.
 

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
	double mean = sum/N, sigma = sqrt(sum2/N - mean*mean);
	//double sigma = sqrt(sum2/N - sum*sum/(N*N));	
	//fprintf(stderr,"sigma is %g",sigma);
	*result = (V) * mean;
	*error = sigma*V/sqrt(N);
}

complex plain(int dim, double f(int dim,double* x), double* a, double* b, int N){
	//The variables
	double V=1, sum=0, sum2=0, x[dim],mean=0,sigma=0;
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

double corput(int n, int base){
	double q=0,bk=(double)1/base;
	while(n>0){q+=(n % base)*bk;n /= base; bk /= base;}
	return q;}

double halton(int n,int k, int d){
	int base[]={2,3,5,7,111,13,17,19,23,29,31,37,41,43,47,53,59,61,67};
	int maxd = sizeof(base)/sizeof(int);
	assert(d<=maxd);
	return corput(n,base[k]);
}

void quasimc(int dim, double f(int dim, double* x),double* a, double* b, int N, double* result, double* error){
//	fprintf(stderr,"HEJ1\n");
	double x[dim], V = 1, sum1 = 0, sum2 =0;
	for(int i =0;i<dim;i++){
		V*=b[i]-a[i];
	}
//	fprintf(stderr,"HEJ\n");
	for(int i=0;i<N/2;i++){
//		fprintf(stderr,"HEJ2\n");
		for(int k=0;k<dim;k++){
	//		fprintf(stderr,"HALTON %g\n",halton(i+1,dim+1,dim+1));
			x[k] = a[k] + halton(i+1,dim+1,dim+1) * (b[k]-a[k]);
//			fprintf(stderr,"%g\n",x[k]);
			}
//		fprintf(stderr,"HEJ3\n");
		sum1+=f(dim,x);
	}
	//fprintf(stderr,"HEJ\n");
	for(int i=0;i<N/2;i++){
		for(int k=0;k<dim;k++){
	//		fprintf(stderr,"HALTON %g\n",halton(i,dim,dim));
			x[k] = a[k] + halton(i,dim,dim) * (b[k]-a[k]);
		}
		sum2+=f(dim,x);
	}
	*result = (sum1+sum2)/N*V;
	*error = fabs(sum1-sum2)/N*V;
//	fprintf(stderr,"result %g and error %g\n",*(double*)result,*(double*)error);
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
	double sum = 1/((M_PI*M_PI*M_PI)*(1 - cos(x[0]) * cos(x[1]) * cos(x[2])));
	return sum;
}

int main(){
	int k = 1;
	double a[k];
	double b[k];
	printf("Part A\n\n");
	printf("The values of the lists that we are integrating x**2 over as a test.\n");
	for(int i=0;i<k;i++){
		a[i] = 0;
		b[i] = 3;
		printf("a[%i]=%g and b[%i]=%g\n",i,a[i],i,b[i]);
	}
	double error=0, result=0;
	int N = 10e6;
	plainmc(k,a,b,f,N,&result,&error);
	printf("The plain integration of x**2 gives %g.\nWith an estimated error of %g\nand an actual error of %g\n",result,error,fabs(result-9));
	// Now we calculate that crazy equation
	int n = 3;
	double a1[] = {0,0,0};
	double b1[] = {M_PI,M_PI,M_PI};
	printf("We now calculate the given integral. \nThe values of the lists that we are integrating over are\n");
	for(int i=0;i<n;i++){
		printf("a[%i]=%g and b[%i]=%g\n",i,a1[i],i,b1[i]);
	}
	double error1=0, result1=0;
	int N1 = 10e5;
	plainmc(n,a1,b1,f1,N1,&result1,&error1);
	printf("montecarlo plain integration of the given function gives %g.\nWith an estimated error of %g\nand an actual error of %g ",result1,error1,fabs(result1-1.3932039296856768591842462603255));
	//Part B
	printf("\n\nPart B\n\n");
	//resetting the error and the result
	error = 0;
	result = 0;
	//We test it by checking if it can solve the integral of x**2 from 0 to 3 in 1 dimension. 
	//This is done because it simply gives 9 and it is thus easy to see if the result is correct.
	for(int i =0; i<k;i++){
		a[i] = 0;
		b[i] = 3;
	}
	quasimc(k,f,a,b,N1,&result,&error);
	printf("We test the quasi-random montecarlo integration of the same function as the one we used to test the pseudo-random integration,\nwhich is the integral of x**2 from 0 to 3.\nFor the quasi-random montegral the integral is %g\nwith an estiamted error of %g\nand an actual error of %g",result,error,fabs(result-9));
	FILE* data = fopen("plot.data.txt","w");
	double actual_q=0;
	double actual_p=0;
	N = 1e6;
//	fprintf("\n\nex 0\n");
	// ----- k=1
	for(int N1=1e5;N1<N;N1+=25000){
		error =0; result = 0; result1 = 0; error1=0;
		plainmc(k,a,b,f,N1,&result1,&error1);
//		fprintf(stderr,"%.10g\n",result1);
		quasimc(k,f,a,b,N1,&result,&error);
//		fprintf(stderr,"%.10g\n",result);
		actual_q = fabs( 9 - result);
		actual_p = fabs( 9 - result1);
		//			n  plainmc  quasimc
		fprintf(data,"%d %.20g %.20g\n",N1,actual_q,actual_p);
		}
	fclose(data);
	return 0;
}

