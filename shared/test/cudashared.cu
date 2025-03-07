#include"gpu.h"

__global__ void neighbor(float *U,float *UR,float *UL);
__global__ void Unew(float *U,float *UR,float *UL);





__device__ float Rusanov(float U_right,float U_left ,float F_right ,float F_left,float speed){

  
return 0.5*(F_right+F_left)-speed*(U_right-U_left);

}

void Call_Unew(){
int threadsPerBlock =512 ;
int blocksPerGrid =(N+threadsPerBlock-1)/threadsPerBlock ;
//size_t size;
Unew<<<blocksPerGrid, threadsPerBlock>>>(d_U,d_UR,d_UL);

}



__global__ void Unew(float *U,float *UR,float *UL){
int I= threadIdx.x;
int i= blockDim.x * blockIdx.x +I ;
float right,left;
float wave;
float FR[2],FL[2];
float Velc,VelR,VelL;
__shared__ float UR1_s[512],UR2_s[512],UL1_s[512],UL2_s[512],U1_s[512],U2_s[512];

	UR1_s[I]=UR[2*i];
        UR2_s[I]=UR[2*i+1];
        U1_s[I]=U[2*i];
        U2_s[I]=U[2*i+1];
        UL1_s[I]=UL[2*i];
        UL2_s[I]=UL[2*i+1];


	if(I<128){
		
		VelR=UR2_s[I]/UR1_s[I];
		Velc=U2_s[I]/U1_s[I];
		VelL=UL2_s[I]/UL1_s[I];
		
	//right
		wave=sqrtf((UR1_s[I]+U1_s[I])*0.5*g);
		right=UR1_s[I]*VelR   ;left=U1_s[I]*Velc     ;
		FR[0]=Rusanov(UR1_s[I],U1_s[I],right,left,wave);
		right=UR1_s[I]*(VelR*VelR+0.5*g*UR1_s[I])   ;left=U1_s[I]*(Velc*Velc+0.5*g*U1_s[I])     ;
		FR[1]=Rusanov(UR2_s[I],U2_s[I],right,left,wave);
	//left
		wave=sqrtf((U1_s[I]+UL1_s[I])*0.5*g);
		right=U1_s[I]*Velc   ;left=UL1_s[I]*VelL     ;
		FL[0]=Rusanov(U1_s[I],UL1_s[I],right,left,wave);
		right=U1_s[I]*(Velc*Velc+0.5*g*U1_s[I])   ;left=UL1_s[I]*(VelL*VelL+0.5*g*UL1_s[I])     ;
		FL[1]=Rusanov(U2_s[I],UL2_s[I],right,left,wave);

	
		U1_s[I]=U1_s[I]-(DT/DX)*(FR[0]-FL[0]);
		U2_s[I]=U2_s[I]-(DT/DX)*(FR[1]-FL[1]);


	}

	U[2*i]=U1_s[I];
	U[2*i+1]=U2_s[I];
}

void Call_neighbor(){
int threadsPerBlock =128 ;
int blocksPerGrid =(N+threadsPerBlock-1)/threadsPerBlock ;
//size_t size;
neighbor<<<blocksPerGrid, threadsPerBlock>>>(d_U,d_UR,d_UL);

}
__global__ void neighbor(float *U,float *UR,float *UL){

int i= blockDim.x * blockIdx.x +threadIdx.x ;
//__shared__ float UR1_s[128],UR2_s[128],UL1_s[128],UL2_s[128],U1_s[128],U2_s[128];


       // UR1_s[I]=UR[2*i];
       // UR2_s[I]=UR[2*i+1];
        //U1_s[I]=U[2*i];
      //  U2_s[I]=U[2*i+1];
    //    UL1_s[I]=UL[2*i];
  //      UL2_s[I]=UL[2*i+1];
  //      __syncthreads();	

	if(i<N){
	

		if(i<(N-1)){
			UR[2*i]=U[2*(i+1)];
			UR[2*i+1]=U[2*(i+1)+1];
	        
		}else{
			UR[2*i]=U[2*i];
			UR[2*i+1]=-U[2*i+1];
		}
		if(i>0){
			UL[2*i]=U[2*(i-1)];
			UL[2*i+1]=U[2*(i-1)+1];
		}else{
			UL[2*i]=U[2*i];
                        UL[2*i+1]=-U[2*i+1];
	
		}

	}

}


void Free_Memory(){

if (h_U) free(h_U);
if (d_U) cudaFree(d_U);
if (d_U) cudaFree(d_UR);
if (d_U) cudaFree(d_UL);
}

void Allocate_Memory(){
size_t size ;
cudaError_t Error ;
size=2*N*sizeof(float);
   h_P=(float*)malloc(size);
   h_U= (float*)malloc(size);

   Error =cudaMalloc((void**)&d_U,size);
   printf("CUDA error (malloc d_U)=%s\n",cudaGetErrorString(Error));
    Error =cudaMalloc((void**)&d_UR,size);
   printf("CUDA error (malloc d_UR)=%s\n",cudaGetErrorString(Error));
 Error =cudaMalloc((void**)&d_UL,size);
   printf("CUDA error (malloc d_UL)=%s\n",cudaGetErrorString(Error));


}

void Send_To_Device(){
size_t size ;
cudaError_t Error ;

size=2*N*sizeof(float);


 Error = cudaMemcpy(d_U,h_U,size,cudaMemcpyHostToDevice);
         printf("CUDA error(memcpy h_U -> d_U)=%s\n",cudaGetErrorString(Error));

}

void Get_From_Device(){

size_t size ;
cudaError_t Error ;


  size=2*N*sizeof(float);


  Error= cudaMemcpy(h_U,d_U,size,cudaMemcpyDeviceToHost);
  printf("CUDA error(memcpy d_U->h_U)=%s\n",cudaGetErrorString(Error));


}
