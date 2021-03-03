#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_rng.h>
#include<pthread.h>
#include<math.h>
//#include<>

struct parameters{int incirle, attempts; const gsl_rng * seed;};

void* my_pi(void* arg){
	struct parameters * data = (struct paramters*)arg;
	(*data).incirle=0;
	double norm, xpoint, ypoint;
	int N=(*data).attempts;
	for(int i=0;i<N;i++){
		xpoint = gsl_rng_uniform(data->seed);
		ypoint = gsl_rng_uniform(data->seed);
		norm = sqrt(pow(xpoint,2) + pow(ypoint,2));
		if(norm<=1){
			(*data).incirle++;
		}
	}
	return NULL;
}

int main(){
return 0;
}
