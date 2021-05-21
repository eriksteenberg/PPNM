#include<gsl/gsl_linalg.h>
#include<stdio.h>
#include<math.h>
#include <gsl/gsl_integration.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_interp.h>
#define RND (double)rand()/RAND_MAX

int binsearch(gsl_vector* x, double z){
	int i=0, j=x->size-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>gsl_vector_get(x,mid)) i=mid; else j=mid;
	}
	return i;
}

double interp(gsl_vector* x, gsl_vector* y,double z){
	int x_len = x->size;
	assert(z>=gsl_vector_get(x,0));
	assert(z<=gsl_vector_get(x,x_len-1));
	assert(x_len >1);	
	int i = binsearch(x,z);
	double slope = (gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/(gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
	double value = gsl_vector_get(y,i) + slope*(z-gsl_vector_get(x,i));
	return value;
}

double interp_quad(gsl_vector* x, gsl_vector* y,double z){
	int x_len = x->size;//gsl_vector_get(x,x_len-1))int x_len = x->size
	assert(z>=gsl_vector_get(x,0));
	assert(z<=gsl_vector_get(x,x_len-1));
	assert(x_len >1);
	int i = binsearch(x,z);
	gsl_vector* c = gsl_vector_alloc(x_len);
	gsl_vector_set(c,1,0);
	double dx(int z){
		return (gsl_vector_get(x,z+1)-gsl_vector_get(x,z));
	}
	double dy(int i){
		return (gsl_vector_get(y,i+1)-gsl_vector_get(y,i));
	}
	double p(int i){
		return dy(i)/dx(i);
	}
	double ci1(gsl_vector* c, int i){
		return (1/dx(i+1)) * (p(i+1) - p(i) - gsl_vector_get(c,i)*dx(i));
	}
	double ci(gsl_vector* c, int i){
		return (1/dx(i)) * (p(i+1) - p(i) - gsl_vector_get(c,i+1)*dx(i+1));
	}
	for(int i=1;i<x_len-2;i++){
		int n = x_len-2;
		gsl_vector_set(c,i+1,ci1(c,i));
	}
	double cin = gsl_vector_get(c,x_len-2);
	gsl_vector_set(c,x_len-2,0.5*cin);
	for(int i=1;i<x_len-2;i++){
		int n = x_len-2;
		gsl_vector_set(c,n-i,ci(c,n-i));
	}
	double b(int i){
		return p(i)-gsl_vector_get(c,i)*dx(i);
	}
	return 	gsl_vector_get(y,i) + b(i)*(z - gsl_vector_get(x,i)) + gsl_vector_get(c,i) * pow(z - gsl_vector_get(x,i),2);
	gsl_vector_free(c);
}
double interp_quad_diff(gsl_vector* x, gsl_vector* y,double z){
	int x_len = x->size;//gsl_vector_get(x,x_len-1))int x_len = x->size
	assert(z>=gsl_vector_get(x,0));
	assert(z<=gsl_vector_get(x,x_len-1));
	assert(x_len >1);
	int i = binsearch(x,z);
	gsl_vector* c = gsl_vector_alloc(x_len);
	gsl_vector_set(c,1,0);
	double dx(int z){
		return (gsl_vector_get(x,z+1)-gsl_vector_get(x,z));
	}
	double dy(int i){
		return (gsl_vector_get(y,i+1)-gsl_vector_get(y,i));
	}
	double p(int i){
		return dy(i)/dx(i);
	}
	double ci1(gsl_vector* c, int i){
		return (1/dx(i+1)) * (p(i+1) - p(i) - gsl_vector_get(c,i)*dx(i));
	}
	double ci(gsl_vector* c, int i){
		return (1/dx(i)) * (p(i+1) - p(i) - gsl_vector_get(c,i+1)*dx(i+1));
	}
	for(int i=1;i<x_len-2;i++){
		int n = x_len-2;
		gsl_vector_set(c,i+1,ci1(c,i));
	}
	double cin = gsl_vector_get(c,x_len-2);
	gsl_vector_set(c,x_len-2,0.5*cin);
	for(int i=1;i<x_len-2;i++){
		int n = x_len-2;
		gsl_vector_set(c,n-i,ci(c,n-i));
	}
	double b(int i){
		return p(i)-gsl_vector_get(c,i)*dx(i);
	}
	return  b(i) + 2* gsl_vector_get(c,i)*(z - gsl_vector_get(x,i));
	gsl_vector_free(c);
}

double interp_quad_integ(gsl_vector* x, gsl_vector* y, double z){
	int x_len = x->size;
	assert(z>=gsl_vector_get(x,0));
	assert(z<=gsl_vector_get(x,x_len-1));
	assert(x_len >1);
	int i = binsearch(x,z);
	gsl_vector* c = gsl_vector_alloc(x_len);
	gsl_vector_set(c,0,0);
	double dx(int z){
		return (gsl_vector_get(x,z+1)-gsl_vector_get(x,z));
	}
	double dy(int i){
		return (gsl_vector_get(y,i+1)-gsl_vector_get(y,i));
	}
	double p(int i){
		return dy(i)/dx(i);
	}
	double ci1(gsl_vector* c, int i) {
		return (1/dx(i+1))*(p(i+1)-p(i)-gsl_vector_get(c,i)*dx(i));
	}
	double ci(gsl_vector* c,int i){
		return (1/dx(i))*(p(i+1)-p(i)-gsl_vector_get(c,i+1)*dx(i+1));
	}
	for (int i=1;i<x_len-2;i++){
		int n = x_len-2;
		gsl_vector_set(c,i+1, ci1(c,i));
	}
	double cin = gsl_vector_get(c, x_len-2);
	gsl_vector_set(c,x_len-2, 0.5*cin);
	for (int i=1;i<x_len-2;i++){
		int n = x_len-2;
		gsl_vector_set(c, n-i,ci(c,n-i));
	}
	double b(int i){
		return p(i)-gsl_vector_get(c,i)*dx(i);
	}
	double ypart = gsl_vector_get(y,i)*(z- gsl_vector_get(x,0));
	double bpart = b(i)*(pow(z,2) /2 - gsl_vector_get(x,i)*z -pow( gsl_vector_get(x,0),2)/2 -  gsl_vector_get(x,i)*gsl_vector_get(x,0));
	double cpart = gsl_vector_get(c,i)*(pow(z,3)/3 - gsl_vector_get(x,i)*pow(z,2) + pow(gsl_vector_get(x,i),2)*z - pow( gsl_vector_get(x,0),3)/3 +  gsl_vector_get(x,i)* pow(gsl_vector_get(x,0),2) - pow( gsl_vector_get(x,i),2)* gsl_vector_get(x,0));
	return ypart + bpart + cpart;
	gsl_vector_free(c);
}
double interp_integ(gsl_vector* x, gsl_vector* y,double z){
	int x_len = x->size;
	assert(z<=gsl_vector_get(x,x_len-1));
	assert(x_len >1);
	int i = binsearch(x,z);			 
	double slope = (gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/(gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
	double integ = gsl_vector_get(y,0)*(z - gsl_vector_get(x,0)) + 0.5*slope*pow((z - gsl_vector_get(x,0)),2);
	return integ;
}
double function(double z){
	double n = 0.5;	
	double resul = n*(z);///(z+1));
	return resul;					 								
}

int main(){
	int n=6;
	gsl_interp* linear = gsl_interp_alloc(gsl_interp_linear,n);
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	double xa[n],ya[n];
	for(int i=0; i< n; i++)
	{
		double j = (double)i;
		gsl_vector_set(x,i,j);
		xa[i]=j;
		double bi=function(j);
		ya[i]=bi;
		gsl_vector_set(y,i,bi);
	}
	double  n1=5;
	int x_len = x->size;
	gsl_vector_get(x,x_len-1);
	gsl_interp_init(linear,xa,ya,n);
	printf("\n\n");
	printf("# index 0: the function\n");
	for(double i=0;i<n1;i+=0.1)
	{
		printf("%10g %10g\n",i,function(i)); 
	}
	printf("\n\n");
	printf("# index 1: Linear. First interpret and then integrate\n");
	for(double i=0;i<n1;i+=0.1)
	{
		double b = interp(x,y,i);
		double integ = interp_integ(x,y,i);
		printf("%10g %10g %10g\n",i,b,integ);
	}
	printf("\n\n");
	printf("# index 2: gsl, first the interpret and then the intergration\n");
	for(double i=0;i<n1;i+=0.1)
	{
		double gsl_l = gsl_interp_eval(linear,xa,ya,i,NULL);
		double gsl_integ = gsl_interp_eval_integ(linear,xa,ya,xa[0],i,NULL);
		printf("%10g %10g %10g\n",i,gsl_l,gsl_integ);
	}
	printf("\n\n");
	printf("# index 3: quad, first interpret and then integration\n");
	for(int i=0; i< n; i++)
	{
		double j = (double)i;
		xa[i]=j;
		double bi=function(j);
		ya[i]=bi;
		gsl_vector_set(y,i,bi);
	}
	for(double i =0;i<n1;i+=0.1)
	{
		double qinterp = interp_quad(x,y,i);
		double qinteg = interp_quad_integ(x,y,i);
		printf("%10g %10g %10g\n",i,qinterp,qinteg);
	}
	printf("\n\n");
	printf("# index 4: quad diff\n");
	for(double i =0;i<n1;i+=0.1)
	{
		double hej = interp_quad_diff(x,y,i);
		printf("%10g %10g\n",i,hej);
	}
	gsl_vector_free(x);
	gsl_vector_free(y);
	return 0;
}
