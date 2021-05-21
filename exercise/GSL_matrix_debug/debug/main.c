#include"stdio.h" 	// to include something you should <> not \"\"
#include"gsl_matrix.h"  // to include something you should <> not \"\" and when calling a gsl header you must start it with gsl/"the header"

int print_half_00(gsl_matrix* m)
{
	double half = 1/2;
	int *status = printf( "half m_{00} = %i\n", gsl_matrix_get(&m,0,0)*half ); //do your calcualtion outside the printf nto inside, by using %i printf prepared for an integer but it receives a double use %g instead. Alo remove the * and & and print doesnt return a int* and &m is a wrong reference
	gsl_matrix_free(m);		//although you can close m inside the funktion I would recommend that you free directly inside main. You must remove one of them.
	return status; //replace with 0, using error number is not recommend for the return command.
}

int main(void)
{
	gsl_matrix m = gsl_matrix_alloc(0,0); // 0,0 creates a null matrix use something else like 1,1 
	gsl_matrix_set(&m,0,0,66); //remove & we are altering m not refering to its address
	printf("half m_{00} (should be 33):\n");
	int *status = print_half_00(&m); // remove * pritn_half returns an int not an int* and the the function is not designed to take the pointer and remove & as we are not printing the address of m
	//It isn't neccesary but I would recommend that you put {} arounf the printf in your if/else commands, just so that you know pricely what is part of the command
	if(status>0)
	printf("status=%g : SOMETHING WENT TERRIBLY WRONG (status>0)\n",*status);
	else
	printf("status=%g : everything went just fine (status=0)\n",*status);
	gsl_matrix_free(&m); //remove & we are remove m not refering to it and decide where you want to free the m, here or in print_half_00
return 0;
}
