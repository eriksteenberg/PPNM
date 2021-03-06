#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g ",gsl_vector_get(v,i));
	printf("\n");
}
/*void matrix_print(char s[], gsl_matrix* v){
        for(int i=0;i<v->size1;i++);
	printf("\n");
	for(int j=0;j<v->size2;j++);
	{printf("%10g ",gsl_matrix_get(v,i,j));
	}
}
*/
//gsl_matrix_get(m,i,j)
int main(){
	int n = 3; //dimension
	gsl_matrix* A = gsl_matrix_alloc(n,n);
//	gsl_matrix*
	gsl_matrix* Acopy = gsl_matrix_alloc(n,n);
	gsl_matrix* Bcopy = gsl_matrix_alloc(n,n);
	gsl_matrix* B = gsl_matrix_alloc(n,n);
	gsl_vector* x = gsl_vector_calloc(n);
	gsl_vector* b = gsl_vector_calloc(n);
//	gsl_vector*i
	gsl_matrix_set_all(A,0);
	gsl_vector* y=gsl_vector_calloc(n);
	gsl_matrix_set_all(Acopy,0);
	//matrix elements
	double aarray[3] = {6.13,-2.90,5.86};
	double barray[3] = {8.08, -6.31, -3.89};
	double carray[3] = {-4.36,1.00,0.19};
	//b vector elements
	double darray[3]={6.23,5.37,2.29};
	for( int i=0; i<A->size1;i++){
		for( int j=0; j<A->size2;j++)
		{
		if(i=0){
		double Aij=aarray[j];
		gsl_matrix_set(A,i,j,Aij);
		}
		if(i=1){
		double Aij=barray[j];
		gsl_matrix_set(A,i,j,Aij);
		}
		if(i=2){
		double Aij=carray[j];
		gsl_matrix_set(A,i,j,Aij);
		}
		}
	}
	for(int i=0; i< b->size; i++)
		{
		double bi=darray[i];
		gsl_vector_set(b,i,bi);
		}
	gsl_matrix_memcpy(Acopy,A);
//	gsl_linalg_HH_solve(Acopy,b,x);
	gsl_matrix_set(B,0,0,aarray[0]);
	gsl_matrix_set(B,0,1,aarray[1]);
	gsl_matrix_set(B,0,2,aarray[2]);
	gsl_matrix_set(B,1,0,barray[0]);
	gsl_matrix_set(B,1,1,barray[1]);
	gsl_matrix_set(B,1,2,barray[2]);
	gsl_matrix_set(B,2,0,carray[0]);
	gsl_matrix_set(B,2,1,carray[1]);
	gsl_matrix_set(B,2,2,carray[2]);
	gsl_matrix_memcpy(Bcopy,B);
	gsl_linalg_HH_solve(Bcopy,b,x);
	vector_print("b is equal to ",b);
	vector_print("Using gsl, I get that x is equal to ",x);
	gsl_blas_dgemv(CblasNoTrans,1,B,x,0,y);
	vector_print("A times x is ",y);
gsl_matrix_free(A);
gsl_matrix_free(B);
gsl_matrix_free(Acopy);
gsl_vector_free(b);
gsl_vector_free(x);
return 0;
}
