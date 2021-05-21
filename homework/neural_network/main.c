#include<stdio.h>
#include<math.h>
#include<stdio.h>
#include"functions.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#define RND (double)rand()/RAND_MAX
double activation_function(double x){
	return x * exp(-x*x);
}

double function(double x){
	return sin(5*x+2)*exp(-x*x);
}

typedef struct{
	int n;
       	double (*f)(double);
       	gsl_vector* params;	
	} ann;

ann* ann_alloc(int n, double (*f)(double)){
	ann* list = malloc(sizeof(ann));
	list->n=n;
	list ->f=f;
	list->params = gsl_vector_alloc(3*n);
	return list;	
}

void ann_free(ann* list){
	gsl_vector_free(list->params);
	free(list);
}
double ann_response(ann* network, double x){
	double sum = 0;
	int N = network->n;
	for(int i=0;i<N;i++){
		double a = gsl_vector_get(network->params,3*i);
		double b = gsl_vector_get(network->params,3*i + 1);
		double w = gsl_vector_get(network->params,3*i + 2);
		sum += network->f((x-a)/b)*w;
	}
	return sum;
}

gsl_vector* input;
gsl_vector* output;
ann* network;

double Cost(gsl_vector* p){
	gsl_vector_memcpy(network->params,p);
	double sum = 0;
	int N = input->size;;
	for(int i =0;i<N;i++){
		double xi = gsl_vector_get(input,i);
		double yi = gsl_vector_get(output,i);
		double fi = ann_response(network,xi);
		sum += fabs(fi-yi)*fabs(fi- yi);
	}
	return sum / N;
}

void ann_train(ann* list, gsl_vector* input, gsl_vector* output){
	int N = list->params->size;
	gsl_vector* p = gsl_vector_alloc(N);
	gsl_vector_memcpy(p,list->params);
	quasinewton(Cost,p,1e-3,1000);
	gsl_vector_memcpy(list->params,p);
	gsl_vector_free(p);
}
	

int main(){	
	int N = 7;
	int M = 3*N;
	input = gsl_vector_alloc(M);
	output = gsl_vector_alloc(M);
	network = ann_alloc(N,activation_function);
	double a0 = -1, aN = 1;

	for(int i = 0; i<M;i++){
		double xi = a0 + (aN-a0)*i/(M-1);
		gsl_vector_set(input,i,xi);
		double fi = function(xi);
		gsl_vector_set(output,i,fi);
		fprintf(stderr,"%g %g \n",xi,fi);
	}
	
	gsl_vector_set(network->params,3*N-1,1);
	//fprintf(stderr,"n = %i %g",network->n,gsl_vector_get(network->params,1));
	for(int I = 0;I<network->n;I++){
		gsl_vector_set(network->params,1 * I ,a0 + (aN -a0)*I/(network->n-1));
		gsl_vector_set(network->params,3 * I + 1,I);
		gsl_vector_set(network->params,3 * I + 2,I);
	}
	ann_train(network,input, output);
	FILE* DATA =  fopen("output.txt","w");
	fprintf(DATA,"# - data\n");
	for(double i =a0;i<aN;i+=1.0/64){ 
		double predicted_output = ann_response(network, i);
		fprintf(DATA,"%g %g \n",i,predicted_output);
	}
	fprintf(DATA,"\n\n # - original\n");
	for(double i =0;i<M;i++){	
		fprintf(DATA,"%g %g\n",gsl_vector_get(input,i),gsl_vector_get(output,i));		
	}
	fclose(DATA);
	ann_free(network);
	gsl_vector_free(input);
	gsl_vector_free(output);
	return 0;
}
