#define EXTERNAL
#include"gpu.h"

FILE *fp1=fopen("results.plt","w");
void plot();
void init();

int main() {
int ttt;
float time ;
struct timeval start,end;
 Allocate_Memory();
 
 init();
 
 Send_To_Device();

gettimeofday(&start,NULL);  
for(ttt=0;ttt<NO_STEPS;ttt++){

Call_neighbor();
Call_Unew();


}
gettimeofday(&end,NULL);   
time=((end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec)/1000000.0);
printf("thread number=%d\tComputation Time = %f sec.\n",128,time);
 Get_From_Device();

 plot();
 
 Free_Memory();



}


void init(){

int i ;
float cx;
 
for(i=0;i<N;i++){

cx=(i+0.5)*DX;

	if(cx<0.5*High){
		h_P[2*i]=10.0;
		h_P[2*i+1]=0.0;

	}else{
   	        h_P[2*i]=1.0;
                h_P[2*i+1]=0.0;
	}

}

for(i=0;i<N;i++){

	h_U[2*i]=h_P[2*i];
	h_U[2*i+1]=h_P[2*i]*h_P[2*i+1];


}


}

void plot(){
int i;
for(i=0;i<N;i++){
	fprintf(fp1,"%g\n",h_U[2*i]);


}

}
