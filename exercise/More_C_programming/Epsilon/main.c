#include<limits.h>
#include<float.h>
#include<stdio.h>

int equal(double a, double b, double tau, double eps);

void name_digit(int i);

int main(){
	printf("-----------Epsilon-----------\npart 1)\n\npart I\n\n");
	printf("max_value=%i\n\n While command used\n",INT_MAX);
	int i=1;
        while(i+1>i) {i++;}
        printf("my max is =%i\n\n for command\n\n",i);
	i = 1;
	for(i=1;i<i+1;i++);
        printf("my max is =%i\n\n do while command\n\n",i);
	i=1;
	do{i++;
	}while(i<i+1);

        printf("my max is =%i\n\npart II\n\n INT_MIN",i);

	printf("min_value=%i\n\n While command used\n",INT_MIN);
	i = 0;
	while(i>i-1){ i+=-1;
	};
        printf("my min is =%i\n\n for command used\n\n",i);
	i = 1;
	for(i=1;i>i-1;i--);
	printf("my min is =%i\n\n do while command used\n\n",i);
	i=1;
	do{
		i--;
	}while(i>i-1);
	
	printf("my min is =%i\n\npart III\n\nMachien epsilon\n",i);
	printf("\n The double epsilon is %g\n\n",DBL_EPSILON);
	double e = 1;
	while(1+e!=1){e/=2;
	} e*=2;
	printf("my double machine epsilon is %g when using a while command\n",e);
	for(e = 1;1+e!=1;e/=2){
	}
	e*=2;
	printf("my double machine epsilon is %g when using a while command\n",e); 
	e = 1;
	do{
		e/=2;
	}while(1+e!=1);
       	e*=2;
	printf("my double machine epsilon is %g when using a while command\n",e);
	printf("\nThe float epsilon is %f\n\n",FLT_EPSILON);
        float e1 = 1;
	while(1+e1!=1){e1/=2;
	} e1*=2;
	printf("my float machine epsilon is %f when using a while command\n",e1);
	for(e1 = 1;1+e1!=1;e1/=2){
	}
	e1*=2;
	printf("my float machine epsilon is %f when using a while command\n",e1);
	e1 = 1;
	do{
		e1/=2;
	}while(1+e1!=1);
	e1*=2;
	printf("my float machine epsilon is %f when using a while command\n",e1);
	printf("\nThe long double epsilon is %Lg\n\n",LDBL_EPSILON);
	long double e2 = 1;
	while(1+e2!=1){e2/=2;
	} e2*=2;
	printf("my long double machine epsilon is %Lg when using a while command\n",e2);
	for(e2 = 1;1+e2!=1;e2/=2){
	}
	e2*=2;
	printf("my long double machine epsilon is %Lg when using a while command\n",e2);
	e2 = 1;
	do{
		e2/=2;
	}while(1+e2!=1);
	e2*=2;
	printf("my long double machine epsilon is %Lg when using a while command\n",e2);
	printf("\n\nPart 2 - up and down sum\n\n");
	float sum_up_float = 0;
	float sum_down_float = 0;
	double sum_up_double = 0;
	double sum_down_double = 0;
	float max = 10;
	max = FLT_MAX/3000000000;
	for(float i1=1;i1<max;i1++){
		fprintf(stderr,"%g\n",i1);
		sum_up_float += 1.0f/i1;
		sum_down_float += 1.0f/(max-i1);
	}
	double DBL_max = 10; 
	double max1 = DBL_MAX/300000000;
	fprintf(stderr,"%g\n",max);
	for(double k=1;k<max1;k++){
		fprintf(stderr,"%g\n",k);
		sum_up_double += 1.0f/k;
		sum_down_double += 1.0f/(max-k);
	}
	printf("\nthe sum up float  is %f\n", sum_up_float);
	printf("\nthe sum down float is %f\n", sum_down_float);
	printf("\nTHe difference comes from the fact the first steps of the down flaot will be so \nsmall that it won't surpase that epsilon limit \nand will not register as changing.\nThe sum up part on the other hand will grow quickly at the beginning but \nat a certain point the additions will also be too small \nfor it to register.\n\n");
	printf("the sum of 1/n is famously known for being divergent, but it will comverge for max as max is a finite point.\n");
	printf("the sum up double is %g\n", sum_up_double);
	printf("the sum down double is %g\n", sum_down_double);
	printf("\nSomething something \n\n");
	printf("Part 3\n\nI have implemented an equal command in equal.c \nand will now test it on two examples\n\n");
	double a = 4;
	double b = 5;
	double b1 = 10;
	double tau = 1;
	double eps = 0.5;
	int ab_diff = equal(a,b,tau,eps);
	printf("equal(%g,%g,%g,%g)=%i\n1 means that it is within the relative equallity\n\n",a,b,tau,eps,ab_diff); //1
	eps = 0.5;
	ab_diff = equal(a,b1,tau,eps);
	printf("equal(%g,%g,%g,%g)=%i\n0 means that it is within the relative equallity\n\n",a,b1,tau,eps,ab_diff); //0
	printf("part 4\nWe test it with whether or not it can recognize 5\nIt returns:\n ");
	name_digit(5);	
	printf("\nIt could. YEASSSS!!!!!\n");
return 0;
}
