#include<stdio.h>
#include<assert.h>
int binary(int n, double* x, double z){
	assert(n>1 && x[0]<=z && z<=x[n-1]);
	int i = 0, j = n-1;
	while(j-i>1){
		int mid = (j-i)/2;
		if(z>x[mid]) i=mid; else j = mid;
		}
	return i;
}

int main(){
	double list[] = {0,1,2,3,5,6,7,8,9,10,11,12,13,14,15};
	double z = 2;
	int n = 15;
	int er_n_i_list = binary(n,list,z);
	printf("%d = %d",n,er_n_i_list);
	return 0;
}
