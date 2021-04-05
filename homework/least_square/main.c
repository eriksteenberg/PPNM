#include"functions.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<gsl/gsl_blas.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX

int main(){
	int m = 2;
	double t[] = {1, 2, 3, 4, 6, 9, 10, 13, 15};
	double y[] = {117,100,88,72,53,29.5,25.2,15.2,11.1};
	int n = sizeof(y)/sizeof(y[0]);
	double dy[n];
	double ylog[n];
	double dlogy[n];
	for(int i=0;i<n;i++){
		dy[i] = 0.05 * y[i];
		dlogy[i] = dy[i]/y[i];
		ylog[i]=log(y[i]);
	}
	double func(int i, double x){
		switch(i){
			case 0: return 1; break;
			case 1: return x; break;
			case 2: return x*x; break;
			case 3: return pow(x,3); break;
			case 4: return pow(x,4); break;
			case 5: return pow(x,5); break;
			default: return NAN;
		}
	}
	gsl_vector* lt = gsl_vector_alloc(n);
	gsl_vector* ly = gsl_vector_alloc(n);
	gsl_vector* ldy = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(lt,i,t[i]);
		gsl_vector_set(ly,i,ylog[i]);
		gsl_vector_set(ldy,i,dlogy[i]);
	}

	gsl_matrix* cov = gsl_matrix_alloc(m,m);
	gsl_vector* c = gsl_vector_alloc(m);
	lsfit(m,func,lt,ly,ldy,c,cov);
	gsl_vector* dfit  = gsl_vector_alloc(m);
	for(int k=0;k<m;k++){
		gsl_vector_set(dfit,k,sqrt(gsl_matrix_get(cov,k,k)));
	}

	double fit(double x){
		double sum =0;
		for(int i=0;i<m;i++){
			sum += gsl_vector_get(c,i)*func(i,x);
	}
		return sum;
	}

	FILE* info=fopen("info.txt","w");
	fprintf(info,"the polynomial is ");
	for(int i=0;i<m;i++){
		fprintf(info,"%g*pow(x,%i) ",gsl_vector_get(c,i),i);
		}
	fprintf(info,"\nThe value of c1 is %g \n",gsl_vector_get(c,1));
	double t0 = -log(2)/gsl_vector_get(c,1);


	fprintf(info,"My estimation of the value of the half-life of ThX, based upon the least square fit, is %g days.\nThe value given on wikipedia is %g days.\n\n",t0,3.6319);
	fprintf(info,"the covariance matrix\n");

	for(int i=0;i<m;i++){
		for(int j=0;j<m;j++){
		fprintf(info,"%3.g ",gsl_matrix_get(cov,i,j));
		}
		fprintf(info,"\n");
	}
	fclose(info);	

	FILE* my_result=fopen("result.txt","w");
	fprintf(my_result,"# data points\n");
	for(int i=0;i<n;i++){
	fprintf(my_result,"%g %g %g\n",t[i],ylog[i],dlogy[i]);
	}
	fprintf(my_result,"\n\n# fit\n");
	for(int i=0;i<n;i++){
	fprintf(my_result,"%g %g \n",t[i],fit(t[i]));
	}
	fprintf(my_result,"\n\n# fit with uncertainties= +delta(c0)+delta(c1)\n");
	for(int i=0;i<n;i++){
	fprintf(my_result,"%g %g \n",t[i],fit(t[i])+gsl_vector_get(dfit,0)*func(0,t[i])+gsl_vector_get(dfit,1)*func(1,t[i]));
	}
	fprintf(my_result,"\n\n# fit with uncertainties= -delta(c0)-delta(c1)\n");
	for(int i=0;i<n;i++){
	fprintf(my_result,"%g %g \n",t[i],fit(t[i])-gsl_vector_get(dfit,0)*func(0,t[i])-gsl_vector_get(dfit,1)*func(1,t[i]));
	}
	fclose(my_result);
	gsl_matrix_free(cov);
	gsl_vector_free(c);
	gsl_vector_free(dfit);
return 0;
}


