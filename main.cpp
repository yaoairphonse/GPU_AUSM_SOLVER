#define EXTERNAL
#include"gpu.h"

void scanbody();
void init();
void plotslice();
void plot3D();
void finalvalue();
void exam();


FILE *fp4=fopen("X37body.txt","r");
FILE *fp1=fopen("slice.plt","w");
FILE *fp2=fopen("3D.plt","w");
FILE *fp3=fopen("examine.plt","w");

int main() {
struct timeval start,end; 
float time ;
int ttt;

	Allocate_Memory();
	scanbody();

	init();
	Send_To_Device();

gettimeofday(&start,NULL);  
	for(ttt=0;ttt<NO_STEPS;ttt++){

		Call_neighbor();

		Call_Unew();


		}
gettimeofday(&end,NULL);  
time=((end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec)/1000000.0);
printf("thread number=%d\tComputation Time = %f sec.\n",512,time);
	Get_From_Device();

	finalvalue();

	plotslice();
	exam();
	plot3D();

       Free_Memory();


return 0;

}
void init(){
int i,j,k ;
int index;
float cx,cy,cz;
index=0;
	for(k=0;k<NZ;k++){
		for(j=0;j<NY;j++){
			for(i=0;i<NX;i++){
			cx=(i+0.5)*DX;					
			cy=(j+0.5)*DY;
			cz=(k+0.5)*DZ;
			// 	
				if(cy<0.02*High){
				h_P[5*index]=0.0001;
				h_P[5*index+1]=0.0;
				h_P[5*index+2]=7000;
				h_P[5*index+3]=0.0;
 				h_P[5*index+4]=213;
				}else{
                                h_P[5*index]=0.0001;
                                h_P[5*index+1]=0.0;
                                h_P[5*index+2]=0.0;
                                h_P[5*index+3]=0.0;
                                h_P[5*index+4]=213;

				}
			
			
			index++;

			}		
		}	

	}
	for(i=0;i<N;i++){
		h_U[5*i]=h_P[5*i];
 		h_U[5*i+1]=h_P[5*i]*h_P[5*i+1];
 		h_U[5*i+2]=h_P[5*i]*h_P[5*i+2];
 		h_U[5*i+3]=h_P[5*i]*h_P[5*i+3];
 		h_U[5*i+4]=h_P[5*i]*(Cv*h_P[5*i+4]+0.5*(h_P[5*i+1]*h_P[5*i+1]+h_P[5*i+2]*h_P[5*i+2]+h_P[5*i+3]*h_P[5*i+3]));

	}

}



void scanbody(){

int i;
int num;
  for(i=0;i<N;i++){
    fscanf(fp4,"%d,",&num);
     if(num<-0.5){
     h_body[i]=0.0;
     	}else{
              h_body[i]=1.0;
     	}

    // very important
    //      if(h_body[i]>0){
    //           total_cells++;
    //                }
    //                  }
  } 
  fclose(fp4);
  
}
void finalvalue(){

int i;
//int x_cell,y_cell,z_cell;

	for(i=0;i<N;i++){

//	        z_cell =(int)i/(NX*NY);
//         	y_cell =(int)(i-z_cell*NY*NX)/NX;
//         	x_cell = i-y_cell*NX-z_cell*NX*NY;


	h_P[5*i]=h_U[5*i];
        h_P[5*i+1]=h_U[5*i+1]/h_U[5*i];
        h_P[5*i+2]=h_U[5*i+2]/h_U[5*i];
        h_P[5*i+3]=h_U[5*i+3]/h_U[5*i];
        h_P[5*i+4]=((h_U[5*i+4]/h_U[5*i])-0.5*(h_P[5*i+1]*h_P[5*i+1]+h_P[5*i+2]*h_P[5*i+2]+h_P[5*i+3]*h_P[5*i+3]))/Cv;

	}
}
void plotslice() {
int i,j,k;
float cx,cy,cz;
int index ;
float CFL;
printf("plotslice...\n");
 fprintf(fp1,"ZONE  i  =  %d   j  = %d     f=point\n",NX,NY); //regardless,it is 100X100X100
index=0;
   for(k=0;k<NZ;k++){
	for(j=0;j<NY;j++){
		for(i=0;i<NX;i++){
		cx =(i+0.5)*DX;
    		cy =(j+0.5)*DY;
		cz =(k+0.5)*DZ;
	        	if(k==0.5*NZ){
			 CFL=DT*(fabs(h_P[5*index+1])/DX+fabs(h_P[5*index+2])/DY+fabs(h_P[5*index+3])/DZ); 			

			 fprintf(fp1,"%f %f %g %g %g %g %g %g\n",cx,cy,CFL,h_P[5*index],h_P[5*index+1],h_P[5*index+2],h_P[5*index+3],h_P[5*index+4]);
				
				}	
            			    index++; 
			
		}
	}
  }
   
 fclose(fp1);
printf("slice is ok...\n");
}
void exam() {
int i,j,k;
float cx,cy,cz;
int index ;
float CFL;
printf("examineslice...\n");
 fprintf(fp3,"ZONE  i  =  %d   j  = %d     f=point\n",NX,NY); //regardless,it is 100X100X100
index=0;
   for(index=0;index<N;index++){
   	        cz =(int)index/(NX*NY);
         	cy =(int)(index-cz*NY*NX)/NX;
         	cx = index-cy*NX-cz*NX*NY;	 

                      if(cz==0.5*NZ-1){
                         CFL=DT*(fabs(h_P[5*index+1])/DX+fabs(h_P[5*index+2])/DY+fabs(h_P[5*index+3])/DZ);

                         fprintf(fp3,"%f %f %g %g %g %g %g %g\n",cx,cy,CFL,h_U[5*index+4],h_UR[5*index+4],h_UL[5*index+4],h_UT[5*index+4],h_UD[5*index+4]);

                                }
                                    

                }
        
  

 fclose(fp3);
printf("exam is ok...\n");
}
void plot3D(){
int i,j,k;
float cx,cy,cz;
int index ;
float CFL;
printf("plot3D...\n");
 fprintf(fp2,"ZONE  i  =  %d   j  = %d  k = %d   f=point\n",NX,NY,NZ);
index=0;
   for(k=0;k<NZ;k++){
        for(j=0;j<NY;j++){
                for(i=0;i<NX;i++){
                cx =(i+0.5)*DX;
                cy =(j+0.5)*DY;
                cz =(k+0.5)*DZ;
		if(h_body[index]>0.5){
		 fprintf(fp2,"%f %f %f %g %g %g %g %g\n",cx,cy,cz,h_P[5*index],h_P[5*index+1],h_P[5*index+2],h_P[5*index+3],28392.0);		
		}else{
		fprintf(fp2,"%f %f %f %g %g %g %g %g\n",cx,cy,cz,h_P[5*index],h_P[5*index+1],h_P[5*index+2],h_P[5*index+3],h_P[5*index+4]);
		}
		index++;		
		}
	}
   }

fclose(fp2);
printf("3D is ok!\n");
}
