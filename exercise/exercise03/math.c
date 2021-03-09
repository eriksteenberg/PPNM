#include<math.h>
#include<stdio.h>
#include<complex.h>
#define pi 3.1415
int main(){
	double x,y;
	complex a,b,c,d,e;
	y = j1(0.5);
	x = tgamma(5);
	a = csqrt(-2);
	b = cexp(I*pi);
	c = exp(I);
	d = pow(I,exp(1));
	e = pow(I,I);
	float x_float = 1.f/9;
	double x_double = 1.9;
	long double x_long_double = 1.L/9;
	printf("------part 1------\n\n");
	printf("gamma(5)=%g\n",x);
	printf("bessel_1(0.5)=%g\n",y);
	printf("csqrt(-2)=%g+I%g\n",creal(a),cimag(a));
	printf("exp(I*PI)=%g+I%g\n",creal(b),cimag(b));
	printf("exp(I)=%g+I%g\n",creal(c),cimag(c));
	printf("I**e=%g+I%g\n",creal(d),cimag(d));
	printf("I**I=%g+I%g\n\n",creal(e),cimag(e));
	printf("------part 2------\n\n");
	printf("float x_float = %f\n",x_float);
	printf("double x_double = %g\n",x_double);
	printf("long double x_double = %Lg\n\n",x_long_double);
return 0;
}
