#include<math.h>
#include<stdio.h>

int main(int argc, char** argv){
        int items;
        double x;
        FILE* input_file = fopen(argv[1],"r");
	FILE* output_file = fopen(argv[2],"w");
	do{
		items = fscanf(input_file,"%lg",&x);
                fprintf(output_file,"x=%g: sin(x)=%g and cos(x)=%g\n",x,sin(x),cos(x));

	}while(items!=EOF);
fclose(input_file);
fclose(output_file);
return 0;
}
