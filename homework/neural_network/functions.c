#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX
static const double DELTA = 2.22045e-16;

void show_matrix(gsl_matrix* A){
	int n = A->size1; 
	int m = A->size2; 
	printf("\n");
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			double x = gsl_matrix_get(A,i,j);
			printf("%9.3g ",x);
		}
		printf("\n");
	}
	printf("\n");
}

double cdot(gsl_vector* A, gsl_vector* B){
	assert(A->size==B->size);	
	double result=0;
	for(int i = 0;i<A->size;i++){
		double b = gsl_vector_get(A,i)*gsl_vector_get(B,i);
		result += b;	
	}
	return result;
}

void show_vector(gsl_vector* V){
	int n = V->size;
	printf("\n");
	for(int i=0;i<n;i++){
		printf("%g\n",gsl_vector_get(V,i));
	}
	printf("\n");
}
void GS_bak(gsl_matrix* R, gsl_vector* x){
	int m=R->size1;
	for(int i=m-1;i>=0;i--){
		double xi = gsl_vector_get(x,i);
		for(int j=i+1;j<m;j++)
			xi-=gsl_matrix_get(R,i,j)*gsl_vector_get(x,j);
		gsl_vector_set(x,i,xi/gsl_matrix_get(R,i,i));
	}
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x);
	GS_bak(R,x);
}

void GS_decomp(gsl_matrix* A,gsl_matrix* R){
	assert(A->size2 == R->size1);
	int N = A->size1;
	int M = A->size2;
	gsl_vector* ai = gsl_vector_alloc(N);
	gsl_vector* aj = gsl_vector_alloc(N);
	for(int i=0;i<M;i++){
                gsl_matrix_get_col(ai,A,i);
		double Rii = gsl_blas_dnrm2(ai);
		gsl_matrix_set(R,i,i,Rii);
		gsl_vector_scale(ai,1/Rii);
		gsl_matrix_set_col(A,i,ai);
		for(int j=i+1;j<M;j++){
	  		gsl_matrix_get_col(aj,A,j);
			double Rij = cdot(ai,aj);
			gsl_blas_daxpy(-Rij,ai,aj);
			gsl_matrix_set(R,i,j,Rij);
			gsl_matrix_set_col(A,j,aj);
		}
	}
	gsl_vector_free(ai);	
	gsl_vector_free(aj);
}

void GS_inv(gsl_matrix* A,gsl_matrix* Inv){
	int n = A->size2;
	assert(n==A->size1);
	gsl_matrix* I = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_matrix* Q = gsl_matrix_alloc(n,n);
	gsl_matrix_memcpy(Q,A);
	gsl_matrix_set_identity(I);
	gsl_vector* ei = gsl_vector_alloc(n);
	gsl_vector* x = gsl_vector_alloc(n);
	GS_decomp(Q,R);
	for(int i = 0;i<n;i++){
		gsl_matrix_get_col(ei,I,i);
		GS_solve(Q,R,ei,x);
		gsl_matrix_set_col(Inv,i,x);
		}
	gsl_matrix_free(I);
	gsl_matrix_free(R);
	gsl_matrix_free(Q);
	gsl_vector_free(ei);
	gsl_vector_free(x);
}


void matrix_multi(gsl_matrix* A, gsl_matrix* B,gsl_matrix* result){
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A,B,0.0,result);
}

void transposed_multi(gsl_matrix* A, gsl_matrix* result){
		gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0,A,A,0.0,result);
}

void lsfit(int m, double f(int i, double x), gsl_vector* x, gsl_vector* y,gsl_vector* dy,gsl_vector* c, gsl_matrix* S){
	int n = y->size;
	assert(y->size==x->size);
   	assert(c->size==m);
	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_matrix* R = gsl_matrix_alloc(m,m);
	gsl_vector* b = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		double xi = gsl_vector_get(x ,i);
		gsl_vector_set(b,i,gsl_vector_get(y,i)/gsl_vector_get(dy,i));
		for(int k=0;k<m;k++){
			gsl_matrix_set(A,i,k,f(k,xi)/gsl_vector_get(dy,i));
		}
	}
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,A,0,S);
	GS_decomp(A,R);	
	GS_solve(A,R,b,c); //vi gemmer ikke i x for ikke at override
	//we also want a covariance matrix	
	GS_inv(S,S);
	gsl_matrix_free(R);
	gsl_matrix_free(A);
	gsl_vector_free(b);
}
//a jacobi taht doesnt return a vector
void jac_special(double f(gsl_vector* x),gsl_vector* x, gsl_matrix* jac){
	int n=x->size;
	double dx = 0.0001;
	double F = f(x);
	for(int i = 0;i < n; i++){
		double xi = gsl_vector_get(x,i);
		gsl_vector_set(x,i,xi + dx);
		double dF = f(x);
		gsl_matrix_set(jac,i,i,(F-dF)/dx);
		gsl_vector_set(x,i,xi);
	}	
}


void jacobi(gsl_vector* x, void f(gsl_vector* xlist, gsl_vector* flist),gsl_matrix* Jac){
	int n = x->size;
	gsl_vector* dflist = gsl_vector_alloc(n);
	gsl_vector* fx = gsl_vector_alloc(n);
	double dx = 0.0001;
	f(x,fx);
	for (int i = 0; i < n; i++){
		double xi = gsl_vector_get(x,i);
		gsl_vector_set(x,i,xi + dx);
		f(x,dflist);
		//gsl_vector_axpby(1,fx,-1,dflist);
		//i TRIED CALCULATING  df-f this way first, but my pc war being annoying at the time
		gsl_vector_sub(dflist,fx);
		gsl_vector_set(x,i,xi);
      		gsl_vector_scale(dflist,1/dx);
		gsl_matrix_set_col(Jac,i,dflist);
	}
	gsl_vector_free(dflist);
	gsl_vector_free(fx);
}
//A combination of all the GS commands to quickly solve the linear equations
void GS_calculate(gsl_matrix* A, gsl_vector* b, gsl_vector* x){
	int n = A->size1;
	gsl_matrix* Q = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_matrix_memcpy(Q,A);
	GS_decomp(Q,R);
	GS_solve(Q,R,b,x);
	//fprintf(stderr,"\nThis should be the x vector\n");
	//show_vector(x);
	gsl_matrix_free(Q);
	gsl_matrix_free(R);
}

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){ 
	//for functions that takes a vector and returns a vector
	int n= x->size;
	gsl_matrix* jac = gsl_matrix_alloc(n,n);
	gsl_vector* dx = gsl_vector_alloc(n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* lambdax = gsl_vector_alloc(n);
	gsl_vector* lfx = gsl_vector_alloc(n);
	f(x,fx);
	int N = 0;
	while (N<1000){
		f(x,fx);
		jacobi(x,f,jac);
		gsl_vector_scale(fx,-1);
		GS_calculate(jac,fx,dx);
		N++;
		double lambda = 2;
		while(1){	
			lambda *=0.5;
			gsl_vector_scale(dx,lambda);
			f(x,fx);
			double norm_fx = gsl_blas_dnrm2(fx);
			gsl_vector_add(x,dx);
			f(x,lfx);//f(x+l*dx)
			double norm_lfx = gsl_blas_dnrm2(lfx);
			if(norm_lfx<(1-0.5*lambda)*norm_fx){//f(x+dx)>c*f(x)
			       //fprintf(stderr,"\n f(x+dx)<const.* f(x)\n");
			       break;}
			if(lambda>0.015625){
		//		fprintf(stderr,"\n Break because of lambda\n");
				break;}	
		}
	}
	gsl_matrix_free(jac);
	gsl_vector_free(dx);
	gsl_vector_free(fx);
	gsl_vector_free(lambdax);
	gsl_vector_free(lfx);
}

void gradient(double f(gsl_vector* xvector),gsl_vector* x,gsl_vector* grad){
	int n = x->size;
	double dx = sqrt(DELTA);
	double fdx = 0;
	double fx = 0;
	double xi = 0;
	double gradi =0;
	for(int i=0;i<n;i++){
		xi = gsl_vector_get(x,i);
		fx = f(x);
		gsl_vector_set(x,i,xi+dx);
		fdx = f(x);
		gradi = (fdx - fx ) / ( dx );
		gsl_vector_set(grad,i,gradi);
		gsl_vector_set(x,i,xi);
		}
}

void quasinewton(double f(gsl_vector* xvector),gsl_vector* x, double eps,int steps){
	int n = x->size;
	gsl_vector* Dx = gsl_vector_alloc(n);
	gsl_vector* grad = gsl_vector_alloc(n);
	gsl_vector* gradxs = gsl_vector_alloc(n);
	gsl_matrix* B = gsl_matrix_alloc(n,n);
	gsl_matrix* I = gsl_matrix_alloc(n,n);
	gsl_vector* xs = gsl_vector_alloc(n);
	gsl_vector* s = gsl_vector_alloc(n);
	gsl_vector* c = gsl_vector_alloc(n);
	gsl_vector* u = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_matrix_set_identity(B);
	gsl_matrix_set_identity(I);
	double fx = f(x);
	double fxs = 0;
	double sTg=0; 
	gradient(f,x,grad);
	int k = 0;
	while(k<steps){
		k++;
		gsl_blas_dgemv(CblasNoTrans, -1, B, grad,0,Dx); // eq(6)
		double lambda = 1;
		double norm_grad = gsl_blas_dnrm2(grad);
		if(norm_grad<DELTA*gsl_blas_dnrm2(x)){
			printf("The gradient is within our acceptance\nit is %g",norm_grad);
			break;
		}
		double norm_Dx = gsl_blas_dnrm2(Dx);
		if(norm_Dx<eps){
			printf("Dx is within our acceptance\n");
			break;
		}
//		printf("time for lambda \n");
		while(1){
//			printf("lambda time\n");
//			gsl_blas_dgemv(CblasNoTrans,lambda,B,Dx,0,s); // eq(8)
			gsl_vector_memcpy(s,Dx);
			gsl_vector_scale(s,lambda);
//			printf("s is equal to \n");
//			show_vector(s);
			if(gsl_blas_dnrm2(s)==0){
				printf("s is zero\n");
				break;
			}
//			printf("lambda is equal to %g k is equal to %i\n",lambda,k);
			gsl_vector_memcpy(xs,x); 
			gsl_vector_add(xs,s); // x+ s
//			printf("xs is");
//			show_vector(xs);
			fxs = f(xs); // f(x+s)
			gsl_blas_ddot(s,grad,&sTg);
//			printf("sTg is %g\n",sTg);
			if(fxs<fx+0.001*sTg){ //checking for eq(9) condition
//				printf("fxs-fx is small\nit is equal to fx=%g and fxs =%g",fx,fxs);
				break;
				}
			if(lambda<0.01){
				gsl_matrix_set_identity(B);
//				printf("lambda is too small\n");
				break;	
			}
			lambda *=0.5;
			}
			if(gsl_blas_dnrm2(Dx)<1e-10){
//				printf("The norm of dx is minimal\n");
				break;
			}
//		printf("lambda time over\n");
		gradient(f,xs,gradxs);
		gsl_vector_memcpy(y,gradxs); 
		gsl_vector_memcpy(u,s); 
		gsl_vector_sub(y,grad); // y = Grad(x+s) - Grad (x)
		gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u);
		double sTy = 0;
		gsl_blas_ddot(s,y,&sTy);
		if(fabs(sTy)>10e-6){
//			printf("B update!\n");
			double gamma = 0;
			double uTy = 0;
			gsl_blas_ddot(u,y,&uTy); // uTy
			gamma = 0.5 *(uTy)*(1/sTy); // sTy
			gsl_vector_memcpy(c,u); 
			gsl_blas_daxpy(-gamma,s,u); // sTy * a
			gsl_vector_scale(u,1/sTy); // a a == u asT + saT
//			show_vector(u);
			gsl_blas_dger(1,u,s,B);
			gsl_blas_dger(1,s,u,B);
//			printf("B is equal to \n");
//			show_matrix(B);
//			gsl_blas_dger(1,c,s,B);
		}
		gsl_vector_memcpy(x,xs);	
		gsl_vector_memcpy(grad,gradxs);
		fx = fxs;
//		printf("fx is equal to %g\n",fx);	

	}
	printf("The number of steps it took to find the extremum is %i\n",k);
	gsl_matrix_free(I);
	gsl_matrix_free(B);
	gsl_vector_free(Dx);
	gsl_vector_free(grad);
	gsl_vector_free(gradxs);
	gsl_vector_free(xs);
	gsl_vector_free(s);
	gsl_vector_free(c);
	gsl_vector_free(u);
	gsl_vector_free(y);

}	
