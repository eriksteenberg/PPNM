#include<math.h>
#include<stdio.h>
#include<stdlib.h>
int main(int argc, char** argv){
	if (argc<2)fprintf(stderr,"tehre were no arguments\n");
	else{
		for(int i=1;i<argc;i++){
			double x = atof(argv[i]);
			double x_sin = sin(x);
			double x_cos = cos(x);
			printf("x=%g: sin(x)=%g and cos(x)=%g\n",x,x_sin,x_cos);
		}
	}
return 0;
}
