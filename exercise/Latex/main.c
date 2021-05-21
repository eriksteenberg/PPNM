#include<stdio.h>
#include<math.h>

double ex(double);
int main(){
	for(double x=-5;x<=5;x+=1./10)
	printf("%g %g %g\n",x,ex(x),exp(x));
return 0;
}
