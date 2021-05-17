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
	c = cexp(I);
	d = cpow(I,exp(1));
	e = cpow(I,I);
	float x_float = 1.f/9;
	double x_double = 1.9;
	long double x_long_double = 1.L/9;
	printf("------part 1------\n\n");
	printf("gamma(5)=%g\n",x);
	printf("bessel_1(0.5)=%g\n",y);
	printf("csqrt(-2)=%g + %gi\n",creal(a),cimag(a));
	printf("exp(i*PI)=%g + %gi\n",creal(b),cimag(b));
	printf("exp(i)=%g + %gi\n",creal(c),cimag(c));
	complex f = cexp(I*pi*exp(1)*0.5);
	printf("i**e=%g + %gi \nwhich according to eulers equation is equal to  \nexp(e*i*pi*0.5)= %g +%gi \n",creal(d),cimag(d),creal(f),cimag(f));
	complex g = cexp(- pi*0.5);
	printf("i**i=%g + %gi\nAccording to eulers equation\nIt should be equal to exp(-pi/2) = %g + %gi\n\n",creal(e),cimag(e),creal(g),cimag(g));
	printf("------part 2------\n\n");
	printf("float x_float = %.25g\n",x_float);
	printf("double x_double = %.25lg\n",x_double);
	printf("long double x_double = %.25Lg\n\n",x_long_double);
return 0;
}
