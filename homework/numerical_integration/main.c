#include<math.h>
#include<stdio.h>
#include<assert.h>
#include<gsl/gsl_integration.h>
#include"functions.h"
static int calls;
double adap_int(double f(double x),double a,double b,double delta, double epsilon){
	double f1 = f(a+(b-a)/6);
	double f2 = f(a+2*(b-a)/6);
	double f3 = f(a+4*(b-a)/6);
	double f4 = f(a+5*(b-a)/6);
	double U = (2*f1+f2+f3+2*f4)/6*(b-a);
	double N = (f1+f2+f3+f4)/4*(b-a);
	double error = fabs(U-N);
	if ( error < delta+epsilon*fabs(U)){return U;}
	else{return adap_int(f,a,(a+b)/2,delta/sqrt(2),epsilon)+adap_int(f,(a+b)/2,b,delta/sqrt(2),epsilon);}
}
static double A,B;
static double F(double f(double), double t){
	return f((A+B)/2 + (A-B)/2*cos(t) ) *sin(t) * (B-A)/2;
	}
double subclenshaw(double f(double),double a, double b, double delta, double epsilon,int number){
	assert(number<1500000);
	double f1 = F(f,(a+(b-a)/6));
	double f2 = F(f,(a+2*(b-a)/6));
	double f3 = F(f,(a+4*(b-a)/6));
	double f4 = F(f,(a+5*(b-a)/6));
	double U = (2*f1+f2+f3+2*f4)/6*(b-a);
	double N = (f1+f2+f3+f4)/4*(b-a);
	double error = fabs(U-N);
	if ( error < delta+epsilon*fabs(U)){
		return U;			
	}
	else{
		return subclenshaw(f,a,(a+b)/2,delta/sqrt(2),epsilon,number+1) + subclenshaw(f,(a+b)/2,b,delta/sqrt(2),epsilon,number+1);
	}
}

double clenshaw(double f(double),double a,double b, double delta, double eps){
	A=a;B=b;a=0;b=M_PI;
	int number = 0;
	return subclenshaw(f,a,b,2*delta,2*eps,number);
}
double func(double t){
	calls++;
	return sqrt(t);
}
double func2(double t){
	calls++;
	return 4*sqrt(1-t*t);
}
double func3(double t){
	calls++;
	return 4*t;
}
double func4(double x){
	calls++;
	return 1/sqrt(x);
}
double func5(double x){
	calls++;
	return log(x)/sqrt(x);
}
double func6(double t){
	calls++;
	return 4*sqrt(1-t*t);
}

int main(){
	double a = 0;
	double b = 1;	
	double eps = 1e-4;
	double tol = 1e-4;
	calls=0;
	double res = adap_int(func,a,b,tol,eps);
	calls=0;
	double res1 = adap_int(func2,a,b,tol,eps);
	printf("The absolute accuracy is %g amd the relative accuracy is %g\n",tol,eps);
	printf("The result of the adaptive integration of sqrt(x) is %.10g\n",res);
	printf("The result of the adaptive integration of 4*sqrt(1-x**2) is %.10g\n",res1);
	calls=0;
	double res2 = clenshaw(func4,a,b,tol,eps);
	printf("\n\n\nThe clenshaw integration of 1/sqrt(x) is %.10g. The number of integrand evaluations is %i\n",res2,calls);
	calls=0;
	double res3 = clenshaw(func5,a,b,tol,eps);
	printf("The clenshaw integration of ln(x)/sqrt(x) is %.10g. The number of integrand evaluations is %i\n",res3,calls);
	calls=0;
	double res4 = clenshaw(func6,a,b,tol,eps);
	printf("The clenshaw integration of 4*sqrt(1-x**2) is %.10g. The number of integrand evaluations is %i\n",res4,calls);
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
	double result,error;
	double alpha = 1.0;
	gsl_function F;
	F.function = &func6;
	F.params = &alpha;
	gsl_integration_qags(&F,0,1,0,1e-6,1000,w,&result,&error);
	printf("The gsl integration of 4*sqrt(1-x**2) gives %.10g",result);
	//printf("The result of the Clenshaw curtis for 1/sqrt(x) is %g\nfor ln(x)/sqrt(x) it is %g \nand for 4sqrt(1-x*x) we get %g\n",res2,res3,res4);
	printf("\nTo compare the integration of 4*sqrt(1-x**2).\nAdaptive integration gives %.10g\nThe clenshaw method gives  %.10g\nThe gsl integration gives  %.10g\n",res1,res4,result);
//	printf("THe gsl integration was give an integration workspace with a value of allocation equal 1000. After the itegration it is %g",w);
	gsl_integration_workspace_free(w);	

return 0;
}	
