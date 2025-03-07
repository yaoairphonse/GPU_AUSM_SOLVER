#include"gpu.h"




__global__ void neighbor(float *U,float *UR,float *UL,float *UT,float *UD,float *UF,float *UB,float *body);
__global__ void Unew(float *U,float *UR,float *UL,float *UT,float *UD,float *UF,float *UB);

__device__ float cor_M(float Vin,float Vout,float Tin,float Tout){
float Min,Mout;
float ain,aout;
float M_plus,M_sub;

ain=sqrtf(GAMA*R*Tin);
aout=sqrtf(GAMA*R*Tout);

Min=Vin/ain;
Mout=Vout/aout;

        if(fabs(Min)>1){
                M_plus=0.5*(Min+fabs(Min));
        }else{
                M_plus=0.25*(Min+1)*(Min+1);
        }
        if(fabs(Mout)>1){
                M_sub=0.5*(Mout-fabs(Mout));
        }else{
                M_sub=-0.25*(Mout-1)*(Mout-1);
        }

                return M_plus+M_sub;
}
__device__ float cor_P(float rhoin,float rhoout,float Vin,float Vout ,float Tin,float Tout){

float Min,Mout;
float ain,aout;
float P_plus,P_sub;
float Pin,Pout;

ain=sqrtf(GAMA*R*Tin);
aout=sqrtf(GAMA*R*Tout);

Min=Vin/ain;
Mout=Vout/aout;

Pin=rhoin*R*Tin;
Pout=rhoout*R*Tout;

        if(fabs(Min)>1){
                P_plus=0.5*((Min+fabs(Min))/Min)*Pin;
        }else{
                P_plus=0.25*(Min+1)*(Min+1)*(2-Min)*Pin;
        }
        if(fabs(Mout)>1){
                P_sub=0.5*((Mout-fabs(Mout))/Mout)*Pout;
        }else{
                P_sub=0.25*(Mout-1)*(Mout-1)*(2+Mout)*Pout;
        }

        return P_plus+P_sub;
}

__device__ float flux(float A,float T){
//A is properties V is velocity T is temperature
float face;
float sonic;

sonic=sqrtf(GAMA*R*T);

face=sonic*A;

return face;

}

__device__ float fluxT(float pE,float T,float rho){
//pE mean rho*E
float face2;
float sonic;

sonic=sqrtf(GAMA*R*T);

face2=sonic*(pE+rho*R*T);

return face2;
}


void Call_Unew(){
int threadsPerBlock =512 ;
int blocksPerGrid =(N+threadsPerBlock-1)/threadsPerBlock ;
//size_t size;
Unew<<<blocksPerGrid, threadsPerBlock>>>(d_U,d_UR,d_UL,d_UT,d_UD,d_UF,d_UB);

cudaDeviceSynchronize();
}

__global__ void Unew(float *U,float *UR,float *UL,float *UT,float *UD,float *UF,float *UB){
int i= blockDim.x * blockIdx.x +threadIdx.x ;
//int x_cell,y_cell,z_cell;
//int five=5;
float Mach;
float press;
float Pc[5],PR[5],PL[5],PT[5],PD[5],PF[5],PB[5]; // properties
float FR[5],FL[5],FT[5],FD[5],FF[5],FB[5];
		if(i<N){
				//z_cell =(int)i/(NX*NY);
				//y_cell =(int)(i-z_cell*NY*NX)/NX;
				//x_cell = i-y_cell*NX-z_cell*NX*NY;				
				Pc[0]=U[5*i];   
				PR[0]=UR[5*i];
				PL[0]=UL[5*i]; 
				PT[0]=UT[5*i]; 
				PD[0]=UD[5*i]; 
				PF[0]=UF[5*i]; 
				PB[0]=UB[5*i];
				//ux
				Pc[1]=U[5*i+1]/U[5*i];	
				PR[1]=UR[5*i+1]/UR[5*i];
				PL[1]=UL[5*i+1]/UL[5*i]; 
				PT[1]=UT[5*i+1]/UT[5*i]; 
				PD[1]=UD[5*i+1]/UD[5*i]; 
				PF[1]=UF[5*i+1]/UF[5*i]; 
				PB[1]=UB[5*i+1]/UB[5*i];	
				//uy
				Pc[2]=U[5*i+2]/U[5*i];
				PR[2]=UR[5*i+2]/UR[5*i];
				PL[2]=UL[5*i+2]/UL[5*i]; 
				PT[2]=UT[5*i+2]/UT[5*i]; 
				PD[2]=UD[5*i+2]/UD[5*i]; 
				PF[2]=UF[5*i+2]/UF[5*i];
				PB[2]=UB[5*i+2]/UB[5*i];	
				//uz
				Pc[3]=U[5*i+3]/U[5*i]; 
				PR[3]=UR[5*i+3]/UR[5*i];
				PL[3]=UL[5*i+3]/UL[5*i]; 
				PT[3]=UT[5*i+3]/UT[5*i]; 
				PD[3]=UD[5*i+3]/UD[5*i]; 
				PF[3]=UF[5*i+3]/UF[5*i]; 
				PB[3]=UB[5*i+3]/UB[5*i];
				//temp...
				Pc[4]=((U[5*i+4]/U[5*i])-0.5*(Pc[1]*Pc[1]+Pc[2]*Pc[2]+Pc[3]*Pc[3]))/Cv; 
				PR[4]=((UR[5*i+4]/UR[5*i])-0.5*(PR[1]*PR[1]+PR[2]*PR[2]+PR[3]*PR[3]))/Cv;
				PL[4]=((UL[5*i+4]/UL[5*i])-0.5*(PL[1]*PL[1]+PL[2]*PL[2]+PL[3]*PL[3]))/Cv; 
				PT[4]=((UT[5*i+4]/UT[5*i])-0.5*(PT[1]*PT[1]+PT[2]*PT[2]+PT[3]*PT[3]))/Cv; 
				PD[4]=((UD[5*i+4]/UD[5*i])-0.5*(PD[1]*PD[1]+PD[2]*PD[2]+PD[3]*PD[3]))/Cv; 
				PF[4]=((UF[5*i+4]/UF[5*i])-0.5*(PF[1]*PF[1]+PF[2]*PF[2]+PF[3]*PF[3]))/Cv; 
				PB[4]=((UB[5*i+4]/UB[5*i])-0.5*(PB[1]*PB[1]+PB[2]*PB[2]+PB[3]*PB[3]))/Cv;
				
				
				//right
		Mach=cor_M(Pc[1],PR[1],Pc[4],PR[4]); press=cor_P(Pc[0],PR[0],Pc[1],PR[1],Pc[4],PR[4]);
		FR[0]=0.5*Mach*(flux(UR[5*i],PR[4])+flux(U[5*i],Pc[4]))-0.5*fabs(Mach)*(flux(UR[5*i],PR[4])-flux(U[5*i],Pc[4]));
		FR[1]=0.5*Mach*(flux(UR[5*i+1],PR[4])+flux(U[5*i+1],Pc[4]))-0.5*fabs(Mach)*(flux(UR[5*i+1],PR[4])-flux(U[5*i+1],Pc[4]))+press;
		FR[2]=0.5*Mach*(flux(UR[5*i+2],PR[4])+flux(U[5*i+2],Pc[4]))-0.5*fabs(Mach)*(flux(UR[5*i+2],PR[4])-flux(U[5*i+2],Pc[4]));
		FR[3]=0.5*Mach*(flux(UR[5*i+3],PR[4])+flux(U[5*i+3],Pc[4]))-0.5*fabs(Mach)*(flux(UR[5*i+3],PR[4])-flux(U[5*i+3],Pc[4]));
		FR[4]=0.5*Mach*(fluxT(UR[5*i+4],PR[4],PR[0])+fluxT(U[5*i+4],Pc[4],Pc[0]))-0.5*fabs(Mach)*(fluxT(UR[5*i+4],PR[4],PR[0])-fluxT(U[5*i+4],Pc[4],Pc[0]));
				//left
                Mach=cor_M(PL[1],Pc[1],PL[4],Pc[4]); press=cor_P(PL[0],Pc[0],PL[1],Pc[1],PL[4],Pc[4]);
                FL[0]=0.5*Mach*(flux(UL[5*i],PL[4])+flux(U[5*i],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i],Pc[4])-flux(UL[5*i],PL[4]));
                FL[1]=0.5*Mach*(flux(UL[5*i+1],PL[4])+flux(U[5*i+1],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i+1],Pc[4])-flux(UL[5*i+1],PL[4]))+press;
                FL[2]=0.5*Mach*(flux(UL[5*i+2],PL[4])+flux(U[5*i+2],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i+2],Pc[4])-flux(UL[5*i+2],PL[4]));
                FL[3]=0.5*Mach*(flux(UL[5*i+3],PL[4])+flux(U[5*i+3],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i+3],Pc[4])-flux(UL[5*i+3],PL[4]));
                FL[4]=0.5*Mach*(fluxT(UL[5*i+4],PL[4],PL[0])+fluxT(U[5*i+4],Pc[4],Pc[0]))-0.5*fabs(Mach)*(fluxT(U[5*i+4],Pc[4],Pc[0])-fluxT(UL[5*i+4],PL[4],PL[0]));
				//top
                Mach=cor_M(Pc[2],PT[2],Pc[4],PT[4]); press=cor_P(Pc[0],PT[0],Pc[2],PT[2],Pc[4],PT[4]);
                FT[0]=0.5*Mach*(flux(UT[5*i],PT[4])+flux(U[5*i],Pc[4]))-0.5*fabs(Mach)*(flux(UT[5*i],PT[4])-flux(U[5*i],Pc[4]));
                FT[1]=0.5*Mach*(flux(UT[5*i+1],PT[4])+flux(U[5*i+1],Pc[4]))-0.5*fabs(Mach)*(flux(UT[5*i+1],PT[4])-flux(U[5*i+1],Pc[4]));
                FT[2]=0.5*Mach*(flux(UT[5*i+2],PT[4])+flux(U[5*i+2],Pc[4]))-0.5*fabs(Mach)*(flux(UT[5*i+2],PT[4])-flux(U[5*i+2],Pc[4]))+press;
                FT[3]=0.5*Mach*(flux(UT[5*i+3],PT[4])+flux(U[5*i+3],Pc[4]))-0.5*fabs(Mach)*(flux(UT[5*i+3],PT[4])-flux(U[5*i+3],Pc[4]));
                FT[4]=0.5*Mach*(fluxT(UT[5*i+4],PT[4],PT[0])+fluxT(U[5*i+4],Pc[4],Pc[0]))-0.5*fabs(Mach)*(fluxT(UT[5*i+4],PT[4],PT[0])-fluxT(U[5*i+4],Pc[4],Pc[0]));
				//DOWN
                Mach=cor_M(PD[2],Pc[2],PD[4],Pc[4]); press=cor_P(PD[0],Pc[0],PD[2],Pc[2],PD[4],Pc[4]);
                FD[0]=0.5*Mach*(flux(UD[5*i],PD[4])+flux(U[5*i],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i],Pc[4])-flux(UD[5*i],PD[4]));
                FD[1]=0.5*Mach*(flux(UD[5*i+1],PD[4])+flux(U[5*i+1],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i+1],Pc[4])-flux(UD[5*i+1],PD[4]));
                FD[2]=0.5*Mach*(flux(UD[5*i+2],PD[4])+flux(U[5*i+2],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i+2],Pc[4])-flux(UD[5*i+2],PD[4]))+press;
                FD[3]=0.5*Mach*(flux(UD[5*i+3],PD[4])+flux(U[5*i+3],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i+3],Pc[4])-flux(UD[5*i+3],PD[4]));
                FD[4]=0.5*Mach*(fluxT(UD[5*i+4],PD[4],PD[0])+fluxT(U[5*i+4],Pc[4],Pc[0]))-0.5*fabs(Mach)*(fluxT(U[5*i+4],Pc[4],Pc[0])-fluxT(UD[5*i+4],PD[4],PD[0]));
				//front
                Mach=cor_M(Pc[3],PF[3],Pc[4],PF[4]); press=cor_P(Pc[0],PF[0],Pc[3],PF[3],Pc[4],PF[4]);
                FF[0]=0.5*Mach*(flux(UF[5*i],PF[4])+flux(U[5*i],Pc[4]))-0.5*fabs(Mach)*(flux(UF[5*i],PF[4])-flux(U[5*i],Pc[4]));
                FF[1]=0.5*Mach*(flux(UF[5*i+1],PF[4])+flux(U[5*i+1],Pc[4]))-0.5*fabs(Mach)*(flux(UF[5*i+1],PF[4])-flux(U[5*i+1],Pc[4]));
                FF[2]=0.5*Mach*(flux(UF[5*i+2],PF[4])+flux(U[5*i+2],Pc[4]))-0.5*fabs(Mach)*(flux(UF[5*i+2],PF[4])-flux(U[5*i+2],Pc[4]));
                FF[3]=0.5*Mach*(flux(UF[5*i+3],PF[4])+flux(U[5*i+3],Pc[4]))-0.5*fabs(Mach)*(flux(UF[5*i+3],PF[4])-flux(U[5*i+3],Pc[4]))+press;
                FF[4]=0.5*Mach*(fluxT(UF[5*i+4],PF[4],PF[0])+fluxT(U[5*i+4],Pc[4],Pc[0]))-0.5*fabs(Mach)*(fluxT(UF[5*i+4],PF[4],PF[0])-fluxT(U[5*i+4],Pc[4],Pc[0]));
				//back	
                Mach=cor_M(PB[3],Pc[3],PB[4],Pc[4]); press=cor_P(PB[0],Pc[0],PB[3],Pc[3],PB[4],Pc[4]);
                FB[0]=0.5*Mach*(flux(UB[5*i],PB[4])+flux(U[5*i],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i],Pc[4])-flux(UB[5*i],PB[4]));
                FB[1]=0.5*Mach*(flux(UB[5*i+1],PB[4])+flux(U[5*i+1],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i+1],Pc[4])-flux(UB[5*i+1],PB[4]));
                FB[2]=0.5*Mach*(flux(UB[5*i+2],PB[4])+flux(U[5*i+2],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i+2],Pc[4])-flux(UB[5*i+2],PB[4]));
                FB[3]=0.5*Mach*(flux(UB[5*i+3],PB[4])+flux(U[5*i+3],Pc[4]))-0.5*fabs(Mach)*(flux(U[5*i+3],Pc[4])-flux(UB[5*i+3],PB[4]))+press;
                FB[4]=0.5*Mach*(fluxT(UB[5*i+4],PB[4],PB[0])+fluxT(U[5*i+4],Pc[4],Pc[0]))-0.5*fabs(Mach)*(fluxT(U[5*i+4],Pc[4],Pc[0])-fluxT(UB[5*i+4],PB[4],PB[0]));
		

				//examine that one time step renew a range with five  (5*i)+0~(5*i+4),no neighbor to interfere.
				//renew
				 U[5*i]=U[5*i]-(DT/DX)*(FR[0]-FL[0])-(DT/DY)*(FT[0]-FD[0])-(DT/DZ)*(FF[0]-FB[0]);
				 U[5*i+1]=U[5*i+1]-(DT/DX)*(FR[1]-FL[1])-(DT/DY)*(FT[1]-FD[1])-(DT/DZ)*(FF[1]-FB[1]);
				 U[5*i+2]=U[5*i+2]-(DT/DX)*(FR[2]-FL[2])-(DT/DY)*(FT[2]-FD[2])-(DT/DZ)*(FF[2]-FB[2]);
				 U[5*i+3]=U[5*i+3]-(DT/DX)*(FR[3]-FL[3])-(DT/DY)*(FT[3]-FD[3])-(DT/DZ)*(FF[3]-FB[3]);
				 U[5*i+4]=U[5*i+4]-(DT/DX)*(FR[4]-FL[4])-(DT/DY)*(FT[4]-FD[4])-(DT/DZ)*(FF[4]-FB[4]);

			        

		
		}
}

void Call_neighbor(){
int threadsPerBlock =512 ;
int blocksPerGrid =(N+threadsPerBlock-1)/threadsPerBlock ;
//size_t size;
neighbor<<<blocksPerGrid, threadsPerBlock>>>(d_U,d_UR,d_UL,d_UT,d_UD,d_UF,d_UB,d_body);

cudaDeviceSynchronize();


}

__global__ void neighbor(float *U,float *UR,float *UL,float *UT,float *UD,float *UF,float *UB,float *body){

int i= blockDim.x * blockIdx.x +threadIdx.x ;
int cx,cy,cz;
int five=5;
	if(i<N){
	        cz =(int)i/(NX*NY);
         	cy =(int)(i-cz*NY*NX)/NX;
         	cx = i-cy*NX-cz*NX*NY;	
			if(cx<(NX-1)){
					if(body[i+1]<0.1){
					//air
						UR[5*i]=U[5*i+five];
						UR[5*i+1]=U[5*i+1+five];
						UR[5*i+2]=U[5*i+2+five];
						UR[5*i+3]=U[5*i+3+five];
						UR[5*i+4]=U[5*i+4+five];
							}else{
							//body
										UR[5*i]=U[5*i];
										UR[5*i+1]=-U[5*i+1];
										UR[5*i+2]=U[5*i+2];
										UR[5*i+3]=U[5*i+3];
										UR[5*i+4]=U[5*i+4];
						}
			}else{
			//right boundary,outlet~
			UR[5*i]=U[5*i];
			UR[5*i+1]=U[5*i+1];
			UR[5*i+2]=U[5*i+2];
			UR[5*i+3]=U[5*i+3];
			UR[5*i+4]=U[5*i+4];
			}
			
			
		 	if(cx>0){
				if(body[i-1]<0.1){
					UL[5*i]=U[5*i-five];
					UL[5*i+1]=U[5*i+1-five];
					UL[5*i+2]=U[5*i+2-five];
					UL[5*i+3]=U[5*i+3-five];
					UL[5*i+4]=U[5*i+4-five];
					}else{
								UL[5*i]=U[5*i];
								UL[5*i+1]=-U[5*i+1];
								UL[5*i+2]=U[5*i+2];
								UL[5*i+3]=U[5*i+3];
								UL[5*i+4]=U[5*i+4];
						}
			}else{
			UL[5*i]=U[5*i];
			UL[5*i+1]=U[5*i+1];
			UL[5*i+2]=U[5*i+2];
			UL[5*i+3]=U[5*i+3];
			UL[5*i+4]=U[5*i+4];				
			}
			
		
			if(cy<NY-1){
			
				if(body[i+NX]<0.1){
				    UT[5*i]=U[5*i+five*NX];
					UT[5*i+1]=U[5*i+1+five*NX];
					UT[5*i+2]=U[5*i+2+five*NX];
					UT[5*i+3]=U[5*i+3+five*NX];
					UT[5*i+4]=U[5*i+4+five*NX];
					}else{
							UT[5*i]=U[5*i];
							UT[5*i+1]=U[5*i+1];
							UT[5*i+2]=-U[5*i+2];
							UT[5*i+3]=U[5*i+3];
							UT[5*i+4]=U[5*i+4];
						}	
			}else{
					UT[5*i]=U[5*i];
					UT[5*i+1]=U[5*i+1];
					UT[5*i+2]=U[5*i+2];
					UT[5*i+3]=U[5*i+3];
					UT[5*i+4]=U[5*i+4];			
			
			
			
			}
		 	if(cy>0){
				if(body[i-NX]<0.1){
					UD[5*i]=U[5*i-five*NX];
					UD[5*i+1]=U[5*i+1-five*NX];
					UD[5*i+2]=U[5*i+2-five*NX];
					UD[5*i+3]=U[5*i+3-five*NX];
					UD[5*i+4]=U[5*i+4-five*NX];
					}else{
					UD[5*i]=U[5*i];
					UD[5*i+1]=U[5*i+1];
					UD[5*i+2]=-U[5*i+2];
					UD[5*i+3]=U[5*i+3];
					UD[5*i+4]=U[5*i+4];
						}
			}else{
			//inlet T=213k rho=1e-4
                                        UD[5*i]=U[5*i];
                                        UD[5*i+1]=U[5*i+1];
                                        UD[5*i+2]=U[5*i+2];
                                        UD[5*i+3]=U[5*i+3];
                                        UD[5*i+4]=U[5*i+4];

			}
			if(cz<NZ-1){
				if(body[i+NX*NY]<0.1){
					UF[5*i]=U[5*i+five*NX*NY];
					UF[5*i+1]=U[5*i+1+five*NX*NY];
					UF[5*i+2]=U[5*i+2+five*NX*NY];
					UF[5*i+3]=U[5*i+3+five*NX*NY];
					UF[5*i+4]=U[5*i+4+five*NX*NY];
					}else{
					UF[5*i]=U[5*i];
					UF[5*i+1]=U[5*i+1];
					UF[5*i+2]=U[5*i+2];
					UF[5*i+3]=-U[5*i+3];
					UF[5*i+4]=U[5*i+4];
					}
			}else{
					UF[5*i]=U[5*i];
					UF[5*i+1]=U[5*i+1];
					UF[5*i+2]=U[5*i+2];
					UF[5*i+3]=U[5*i+3];
					UF[5*i+4]=U[5*i+4];
			
			
			}
			if(cz>0){
			
				if(body[i-NX*NY]<0.1){
					UB[5*i]=U[5*i-five*NX*NY];
					UB[5*i+1]=U[5*i+1-five*NX*NY];
					UB[5*i+2]=U[5*i+2-five*NX*NY];
					UB[5*i+3]=U[5*i+3-five*NX*NY];
					UB[5*i+4]=U[5*i+4-five*NX*NY];
					}else{
					UB[5*i]=U[5*i];
					UB[5*i+1]=U[5*i+1];
					UB[5*i+2]=U[5*i+2];
					UB[5*i+3]=-U[5*i+3];
					UB[5*i+4]=U[5*i+4];
						}	
			}else{
					UB[5*i]=U[5*i];
					UB[5*i+1]=U[5*i+1];
					UB[5*i+2]=U[5*i+2];
					UB[5*i+3]=U[5*i+3];
					UB[5*i+4]=U[5*i+4];		
					
			}
	}

}





void Send_To_Device() {
// Size of data to send

size_t size ;
cudaError_t Error ;

	 size = N*sizeof(float);
         Error = cudaMemcpy(d_body,h_body,size,cudaMemcpyHostToDevice);
         printf("CUDA error(memcpy h_body -> d_body)=%s\n",cudaGetErrorString(Error));


	 size = 5*N*sizeof(float);
//	 Error = cudaMemcpy(d_P,h_P,size,cudaMemcpyHostToDevice);
//	 printf("CUDA error(memcpy h_P -> d_P)=%s\n",cudaGetErrorString(Error));
	 Error = cudaMemcpy(d_U,h_U,size,cudaMemcpyHostToDevice);
         printf("CUDA error(memcpy h_U -> d_U)=%s\n",cudaGetErrorString(Error));

}
void Get_From_Device(){


size_t size ;
cudaError_t Error ;

     size = N*sizeof(float);
  Error= cudaMemcpy(h_body,d_body,size,cudaMemcpyDeviceToHost);
  printf("CUDA error(memcpy d_body->h_body)=%s\n",cudaGetErrorString(Error));



  size=5*N*sizeof(float);
  
  Error= cudaMemcpy(h_U,d_U,size,cudaMemcpyDeviceToHost);
  printf("CUDA error(memcpy d_U->h_U)=%s\n",cudaGetErrorString(Error));
//  Error= cudaMemcpy(h_P,d_P,size,cudaMemcpyDeviceToHost);
//  printf("CUDA error(memcpy d_P->h_P)=%s\n",cudaGetErrorString(Error));

    Error= cudaMemcpy(h_UR,d_UR,size,cudaMemcpyDeviceToHost);
  printf("CUDA error(memcpy d_U->h_U)=%s\n",cudaGetErrorString(Error));
  Error= cudaMemcpy(h_UL,d_UL,size,cudaMemcpyDeviceToHost);
  printf("CUDA error(memcpy d_U->h_U)=%s\n",cudaGetErrorString(Error));
  Error= cudaMemcpy(h_UT,d_UT,size,cudaMemcpyDeviceToHost);
  printf("CUDA error(memcpy d_U->h_U)=%s\n",cudaGetErrorString(Error));
  Error= cudaMemcpy(h_UD,d_UD,size,cudaMemcpyDeviceToHost);
  printf("CUDA error(memcpy d_U->h_U)=%s\n",cudaGetErrorString(Error));
  Error= cudaMemcpy(h_UF,d_UF,size,cudaMemcpyDeviceToHost);
  printf("CUDA error(memcpy d_U->h_U)=%s\n",cudaGetErrorString(Error));
  Error= cudaMemcpy(h_UB,d_UB,size,cudaMemcpyDeviceToHost);
  printf("CUDA error(memcpy d_U->h_U)=%s\n",cudaGetErrorString(Error));


}

void Allocate_Memory(){
 
   size_t size ;
   cudaError_t Error ;

   size=N*sizeof(float) ;
     h_body = (float*)malloc(size) ;
    Error =cudaMalloc((void**)&d_body,size);
   printf("CUDA error (malloc d_body)=%s\n",cudaGetErrorString(Error));

   size= 5*N*sizeof(float) ;
   
   

   h_P = (float*)malloc(size) ;  
   h_U = (float*)malloc(size) ;
   h_UR = (float*)malloc(size) ;
   h_UL = (float*)malloc(size) ;
   h_UT = (float*)malloc(size) ;
   h_UD = (float*)malloc(size) ;
   h_UF = (float*)malloc(size) ;
   h_UB = (float*)malloc(size) ;



  // Error =cudaMalloc((void**)&d_P,size);
  // printf("CUDA error (malloc d_P)=%s\n",cudaGetErrorString(Error));
   Error =cudaMalloc((void**)&d_U,size);
   printf("CUDA error (malloc d_U)=%s\n",cudaGetErrorString(Error));
   Error =cudaMalloc((void**)&d_UR,size);
   printf("CUDA error (malloc d_UR)=%s\n",cudaGetErrorString(Error));
   Error =cudaMalloc((void**)&d_UL,size);
   printf("CUDA error (malloc d_UL)=%s\n",cudaGetErrorString(Error));
   Error =cudaMalloc((void**)&d_UT,size);
   printf("CUDA error (malloc d_UT)=%s\n",cudaGetErrorString(Error));
   Error =cudaMalloc((void**)&d_UD,size);
   printf("CUDA error (malloc d_UD)=%s\n",cudaGetErrorString(Error));
   Error =cudaMalloc((void**)&d_UF,size);
   printf("CUDA error (malloc d_UF)=%s\n",cudaGetErrorString(Error));
   Error =cudaMalloc((void**)&d_UB,size);
   printf("CUDA error (malloc d_UB)=%s\n",cudaGetErrorString(Error));


}

void Free_Memory(){
	if (h_body) free(h_body);
	if (d_body) cudaFree(d_body);
        if (h_U) free(h_U);
        if (h_P) free(h_P);
      //  if (d_P) cudaFree(d_P);
	if (d_U) cudaFree(d_U);
	if (d_UR) cudaFree(d_UR);
 	if (d_UL) cudaFree(d_UL);
 	if (d_UT) cudaFree(d_UT);
 	if (d_UD) cudaFree(d_UD);
 	if (d_UF) cudaFree(d_UF);
 	if (d_UB) cudaFree(d_UB);

        if (h_UR) free(h_UR);
        if (h_UL) free(h_UL);
        if (h_UT) free(h_UT);
        if (h_UD) free(h_UD);
        if (h_UF) free(h_UF);
        if (h_UB) free(h_UB);



}
