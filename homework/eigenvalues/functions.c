#include<gsl/gsl_vector.h>
//#include<gsl/gsl_vector_double.h>
//#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
//#include<assert.h>
//#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX

//I've removed and added functions to this functions.c folder, My goal to to create one great functions folder, bu sometimes I remove functions that Arent relevant In a desperate attempt to make my code work

double norm(double x){
	double result = sqrt(pow(x,2));
	return result;
}

void show_matrix(gsl_matrix* A){
	int n = A->size1; 
	int m = A->size2; 
	printf("\n");
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			double x = gsl_matrix_get(A,i,j);
				if(norm(x)<10e-20){	// I added this part because seeing 10-16 5 times, hurts my eyes
					printf("%9i ",0);
				}
				else{
					printf("%9.3g ",x);
				}
		}
		printf("\n");
	}
	printf("\n");
}

void show_vector(gsl_vector* V){
	int n = V->size;
	printf("\n");
	for(int i=0;i<n;i++){
		printf("%g\n",gsl_vector_get(V,i));
	}
	printf("\n");
}
void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c= cos(theta), s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip = -s*gsl_matrix_get(A,i,q) + c*gsl_matrix_get(A,i,p);
		double new_aiq = c*gsl_matrix_get(A,i,q) + s*gsl_matrix_get(A,i,p);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}
		
void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta), s=sin(theta);
	for(int j=0;j<A->size2;j++){
		double new_apj=c*gsl_matrix_get(A,p,j) + s*gsl_matrix_get(A,q,j);
		double new_aqj=-s*gsl_matrix_get(A,p,j) + c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_apj);
		gsl_matrix_set(A,q,j,new_aqj);
	}
}

void jacobi_diag(gsl_matrix* A, gsl_matrix* V){
	int n=A->size1;
	int changed;
	do{
		changed = 0;
		for(int p=0;p<n-1;p++){
			for(int q=p+1;q<n;q++){
				double app=gsl_matrix_get(A,p,p);
				double apq=gsl_matrix_get(A,p,q);
				double aqq=gsl_matrix_get(A,q,q);
				double theta=0.5*atan2(2*apq,(aqq-app));
				double c=cos(theta),s=sin(theta);
				double new_aqq=s*s*app-2*s*c*apq+c*c*aqq;
				double new_app=c*c*app-2*s*c*apq+s*s*aqq;
				if(new_app!=app || new_aqq!=aqq){
					changed=1;
					timesJ(A,p,q,theta);
					Jtimes(A,p,q,-theta);
					timesJ(V,p,q,theta);
				}
			}
		}
	}while(changed!=0);
}
