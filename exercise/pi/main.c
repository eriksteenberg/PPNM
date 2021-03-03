#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<math.h>
#include<omp.h>

struct params {int N, Ninc; unsigned int seed;};

void* my_pi(void* args){
	struct params * data = (struct params*)args;
	(*data).Ninc=0;
	int N=(*data).N;
	for(int i=0;i<N;i++){
		double xpoint = (double) rand_r(&(data->seed))/RAND_MAX;
		double ypoint = (double) rand_r(&(data->seed))/RAND_MAX;
		double norm = sqrt(pow(xpoint,2) + pow(ypoint,2));
		if(norm<=1)(*data).Ninc++;
	}
return NULL;
}

int main(){
	int Ntot = (int)1e8;
	struct params n1 = {.N=Ntot/3, .Ninc=0, .seed = 5};
	struct params n2 = {.N=Ntot/3, .Ninc=0, .seed = 7};
	struct params n3 = {.N=Ntot/3, .Ninc=0, .seed = 15};
	pthread_t p1,p2,p3;
	pthread_create(&p1,NULL,my_pi,(void*)&n1);
	pthread_create(&p2,NULL,my_pi,(void*)&n2);
	pthread_create(&p3,NULL,my_pi,(void*)&n3);
	pthread_join(p1,NULL);
	pthread_join(p2,NULL);
	pthread_join(p3,NULL);
	int Ninc = n1.Ninc + n2.Ninc + n3.Ninc;
	double pi = 4 * (double)Ninc/Ntot;
//	printf("%i %i %10g\n",Ninc,Ntot,(double)Ninc/Ntot);
	printf("my monte carlo simulation of pi gives %g. \nThe value of pi from math.h is %g \n",pi,M_PI);
return 0;
}
