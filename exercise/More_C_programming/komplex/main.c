#include"komplex.h"
#include<stdio.h>

int main(){
	komplex z = {1,1};
	komplex y = {1,-1};
	komplex x = komplex_add(z,y);
//	komplex_set(*z,1,1);
	komplex_print("This is the z",z);
	komplex_print("This is the y",y);
	komplex_print("This is the x, which is the sum of z and y",x);
	komplex b = komplex_new(3,2);
	komplex_print("This is the b which is set to",b);
	komplex d = komplex_sub(y,b);
	komplex_print("this is d, which is y - b",d);
	printf("\nThis complets the illustration of the komplex structure\n");
	return 0;
}
