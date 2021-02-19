#include<math.h>
#include<stdio.h>

int main(){
        int items;
	double x;			        
	do{
		items = fscanf(stdin,"%lg",&x);
		printf("x=%g: sin(x)=%g and cos(x)=%g\n",x,sin(x),cos(x));
	}while(items!=EOF);
return 0;
}



