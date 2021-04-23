#include<math.h>
#include<stdio.h>


double adap_int(double f(double x),double a,double b,double epsilon,double delta){
	double f1 = f(a+(b-a)/6);
	double f2 = f(a+2*(b-a)/6);
	double f3 = f(a+4*(b-a)/6);
	double f4 = f(a+5*(b-a)/6);
	double U = (2*f1+f2+f3+2*f4)/6*(b-a);
	double N = (f1+f2+f3+f4)/4*(b-a);
	double error = fabs(U-N);
	if ( error < delta+epsilon*fabs(U)){
		return U;
	}
	else{
		return adap_int(f,a,(a+b)/2,delta/sqrt(2),epsilon)+ 
			adap_int(f,(a+b)/2,b,delta/sqrt(2),epsilon);
	}
}
double func(double t){
	return sqrt(t);
}
double func2(double t){
	return 4*sqrt(1-t*t);
}
double func3(double t){
	return 4*t;
}
static double A,B;

double g(double f(double x),double x){
	double y = (A + B) / 2 + ((B - A) / 2);
	return f(y * cos(x)) * sin(x) * (B - A) / 2;	
	}

double clenshaw_curtis(double f(double x),double a, double b, double tol, double eps){
	A=a;
	B=b;
	double f1 = g(f,a+(b-a)/6);
	double f2 = g(f,a+2*(b-a)/6);
	double f3 = g(f,a+4*(b-a)/6);
	double f4 = g(f,a+5*(b-a)/6);
	double U = (2*f1+f2+f3+2*f4)/6*(b-a);
	double N = (f1+f2+f3+f4)/4*(b-a);
	double error = fabs(U-N);
	if ( error < tol+eps*fabs(U)){
		return U;
	}
	else{
		return clenshaw_curtis(f,a,(a+b)/2,tol/sqrt(2),eps)+
			clenshaw_curtis(f,(a+b)/2,b,tol/sqrt(2),eps);	

	return adap_int(f,a,b,tol,eps);
	}
}
int main(){
	double a = 0;
	double b = 1;	
	double eps = 1e-6;
	double tol = 1e-6;
	double res = adap_int(func,a,b,tol,eps);
	double res1 = adap_int(func2,a,b,tol,eps);
	fprintf(stderr,"The result of the adaptive integration of sqrt(x) is %g\n and for 4*sqrt(1-x**2) it is %g\n",res,res1);
/*	double res2 = clenshaw_curtis(func,a,b,tol,eps);
	double res3 = 1;//clenshaw_curtis(func2,a,b,tol,eps);
	fprintf(stderr,"The result of the Clenshaw curtis for func1 is %g\n and for func2 is %g \n",res2,res3);
*/return 0;
}	
