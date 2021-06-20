#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX

double norm(double x){
	double result = sqrt(pow(x,2));
	return result;
}
void show_matrix(gsl_matrix* A){
	int n = A->size1; 
	int m = A->size2; 
	fprintf(stderr,"\n");
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			double x = gsl_matrix_get(A,i,j);
				if(norm(x)<10e-10){	// I added this part because seeing 10-16 5 times, hurts my eyes
					fprintf(stderr,"%9i ",0);
				}
				else{
					fprintf(stderr,"%9.3g ",x);
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
		fprintf(stderr,"%g\n",gsl_vector_get(V,i));}
	printf("\n");
}

void rkstep12(void f(double t, gsl_vector* y, gsl_vector* dydt), double t0, gsl_vector* y0, double h, gsl_vector* yend, gsl_vector* yerr){
	int n = y0->size;	//k0 is dydt
	gsl_vector* k0=gsl_vector_alloc(n);	//we generate dydt
	gsl_vector* k1=gsl_vector_alloc(n);	//we generate dydt
	gsl_vector* k2=gsl_vector_alloc(n);	//we generate dydt
	f(t0,y0,k0);	//Now we calculate yend, by euler
	gsl_vector_memcpy(k1,y0);
	gsl_blas_daxpy(0.5*h,k0,k1);
	f(t0+0.5*h,k1,k2);
	for(int i=0;i<n;i++){
		double k1i = gsl_vector_get(k2,i);
		double y0i = gsl_vector_get(y0,i);
		gsl_vector_set(yend,i,y0i+h*k1i);
	}
	double b1_star = 1;
	double b2 = 1;
	for(int i = 0;i<n;i++){
		gsl_vector_set(yerr,i,h*(b2*gsl_vector_get(k2,i)-b1_star*gsl_vector_get(k0,i)));
		}
	gsl_vector_free(k0);
	gsl_vector_free(k1);
	gsl_vector_free(k2);
}

void driver(void f(double t, gsl_vector* y, gsl_vector* dydt),double a,
	double b,double h, gsl_vector* y,gsl_vector* yerr, double acc,double eps){
	FILE* RESULTS = fopen("result.txt","a");
	int n = y->size;
	double err, normy,tol,k=a;
	gsl_vector* yb = gsl_vector_alloc(n);
	while(k<b){
		if(k+h>b)h=b-k;
       		rkstep12(f,k,y,h,yb,yerr);
		err = gsl_blas_dnrm2(yerr);
		normy = gsl_blas_dnrm2(yb);
		tol = (normy * eps + acc) * sqrt(h/(b-a));
		if(err<tol){
			k += h;	
			gsl_vector_memcpy(y,yb);
			//fprintf(RESULTS,"%11g ",k);
			for(int i=0;i<n;i++){
				fprintf(RESULTS,"%11g ",(gsl_vector_get(y,i)));
			}
			fprintf(RESULTS,"\n");

		}
		if(err>0)h*= pow(tol/err,0.25)*0.95;
		else h*=2;
	}
	fprintf(RESULTS,"\n\n");
	fclose(RESULTS);
	gsl_vector_free(yb);
	}

void fun(double t, gsl_vector* y, gsl_vector* dydt){
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
}
void GR(double t, gsl_vector* y, gsl_vector* dydt){
	double m1 = 1;	
	double m2 = 1;	
	double m3 = 1;
	double G = 1;	
	// get s
	double s1x = gsl_vector_get(y,0);
	double s1y = gsl_vector_get(y,1);
	double s2x = gsl_vector_get(y,2);
	double s2y = gsl_vector_get(y,3);
	double s3x = gsl_vector_get(y,4);
	double s3y = gsl_vector_get(y,5);
	//get v
	double v1x = gsl_vector_get(y,6);
	double v1y = gsl_vector_get(y,7);
	double v2x = gsl_vector_get(y,8);
	double v2y = gsl_vector_get(y,9);
	double v3x = gsl_vector_get(y,10);
	double v3y = gsl_vector_get(y,11);
	// x and y distances between bodies
	//fprintf(stderr,"The velocities are equaul to %g %g %g %g %g %g\n",v1x,v1y,v2x,v2y,v3x,v3y);
	double r12x = s2x - s1x; //r21 = -r12
	double r12y = s2y - s1y;
	double r13x = s3x - s1x;
	double r13y = s3y - s1y;
	double r23x = s3x - s2x;
	double r23y = s3y - s2y;
	// set up distances
	double r12 = sqrt(pow(r12x,2)+pow(r12y,2));
	double r13 = sqrt(pow(r13x,2)+pow(r13y,2));
	double r23 = sqrt(pow(r23x,2)+pow(r23y,2));
	//setting up the x and y components of the forces between the bodies
	double f12 = G*m1*m2 /(pow(r12,3));
	double f13 = G*m1*m3 /(pow(r13,3));
	double f32 = G*m2*m3 /(pow(r23,3));
	double f1x = (1/m1)*(f12*(r12x) + f13*r13x);
	double f1y = (1/m1)*(f12*(r12y) + f13*r13y);
	double f2x = (1/m2)*(f12*(-r12x) + f32*r23x);
	double f2y = (1/m2)*(f12*(-r12y) + f32*r23y);
	double f3x = (1/m3)*(f32*(-r23x) + f13*(-r13x));
	double f3y = (1/m3)*(f32*(-r23y) + f13*(-r13y));
	//fprintf(stderr," The forces are equal to %g %g %g %g %g %g\n",f1x,f1y,f2x,f2y,f3x,f3y);
	gsl_vector_set(dydt,0,v1x);
	gsl_vector_set(dydt,1,v1y);
	gsl_vector_set(dydt,2,v2x);
	gsl_vector_set(dydt,3,v2y);
	gsl_vector_set(dydt,4,v3x);
	gsl_vector_set(dydt,5,v3y);
	gsl_vector_set(dydt,6,f1x);
	gsl_vector_set(dydt,7,f1y);
	gsl_vector_set(dydt,8,f2x);
	gsl_vector_set(dydt,9,f2y);
	gsl_vector_set(dydt,10,f3x);
	gsl_vector_set(dydt,11,f3y);
	//fprintf(stderr,"THe y vector\n");
	//show_vector(y);
	}
int main(){
	int n = 13;	// s1 s2 v1 v2 m1 m2 G
	gsl_vector* y0   = gsl_vector_alloc(n);
	gsl_vector* yend = gsl_vector_alloc(n);
	gsl_vector* dydt = gsl_vector_alloc(n);
	gsl_vector* yerr = gsl_vector_alloc(n);
	double dt = 0.1;
	double t0 = 0;
	double tend =6;
	double s1x = -0.97000436;
	double s1y = 0.24308753;
	double s2x = 0;
	double s2y = 0;
	double s3x = 0.97000436;
	double s3y = -0.24308753;
	double v1x = 0.4662036850;
	double v1y = 0.4323657300;
	double v2x = -0.93240737;
	double v2y = -0.86473146;
	double v3x = 0.4662036850;
	double v3y = 0.4323657300;
	// ekstra for dydt
	double F = 0;
	//set up y0
	//set up s
	gsl_vector_set(y0,0,s1x);
	gsl_vector_set(y0,1,s1y);
	gsl_vector_set(y0,2,s2x);
	gsl_vector_set(y0,3,s2y);
	gsl_vector_set(y0,4,s3x);
	gsl_vector_set(y0,5,s3y);
	// set up v
	gsl_vector_set(y0,6,v1x);
	gsl_vector_set(y0,7,v1y);
	gsl_vector_set(y0,8,v2x);
	gsl_vector_set(y0,9,v2y);
	gsl_vector_set(y0,10,v3x);
	gsl_vector_set(y0,11,v3y);
	//set up dydt
	gsl_vector_set(dydt,0,v1x);
	gsl_vector_set(dydt,1,v1y);
	gsl_vector_set(dydt,2,v2x);
	gsl_vector_set(dydt,3,v2y);
	gsl_vector_set(dydt,4,v3x);
	gsl_vector_set(dydt,5,v3y);
	gsl_vector_set(dydt,6,F);
	gsl_vector_set(dydt,7,F);
	gsl_vector_set(dydt,8,F);
	gsl_vector_set(dydt,9,F);
	gsl_vector_set(dydt,10,F);
	gsl_vector_set(dydt,11,F);
	double acc = 1e-2;
	double eps = 1e-2;
	fprintf(stderr,"\n\n");
	FILE* RESULTS = fopen("result.txt","w");
	fprintf(RESULTS,"\n # index 0 : N-body system\n");
	fclose(RESULTS); 
	driver(GR, t0, tend,dt,y0,yerr,acc,eps);
	fprintf(stderr,"\n\n");
	gsl_vector_free(y0);
	gsl_vector_free(yend);
	gsl_vector_free(dydt);
	gsl_vector_free(yerr);
	return 0;
}

