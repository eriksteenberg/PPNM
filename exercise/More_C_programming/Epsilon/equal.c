#include<math.h>
int equal(double a, double b, double tau, double epsilon){
	double abs = fabs(a-b);
	if(abs<tau){
		return 1;
	}
	else if(abs/(fabs(a)+fabs(b))<epsilon*0.5){
		return 1;
	}
return 0;
}
