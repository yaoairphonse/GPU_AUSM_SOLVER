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
int i= blockDim.x * blockIdx.x +threadIdx.x ;
float right,left;
float wave;
float FR[2],FL[2];
float Velc,VelR,VelL;
	if(i<N){
		
		VelR=UR[2*i+1]/UR[2*i];
		Velc=U[2*i+1]/U[2*i];
		VelL=UL[2*i+1]/UL[2*i];
		
	//right
		wave=sqrtf((UR[2*i]+U[2*i])*0.5*g);
		right=UR[2*i]*VelR   ;left=U[2*i]*Velc     ;
		FR[0]=Rusanov(UR[2*i],U[2*i],right,left,wave);
		right=UR[2*i]*(VelR*VelR+0.5*g*UR[2*i])   ;left=U[2*i]*(Velc*Velc+0.5*g*U[2*i])     ;
		FR[1]=Rusanov(UR[2*i+1],U[2*i+1],right,left,wave);
	//left
		wave=sqrtf((U[2*i]+UL[2*i])*0.5*g);
		right=U[2*i]*Velc   ;left=UL[2*i]*VelL     ;
		FL[0]=Rusanov(U[2*i],UL[2*i],right,left,wave);
		right=U[2*i]*(Velc*Velc+0.5*g*U[2*i])   ;left=UL[2*i]*(VelL*VelL+0.5*g*UL[2*i])     ;
		FL[1]=Rusanov(U[2*i+1],UL[2*i+1],right,left,wave);

	
		U[2*i]=U[2*i]-(DT/DX)*(FR[0]-FL[0]);
		U[2*i+1]=U[2*i+1]-(DT/DX)*(FR[1]-FL[1]);


	}
}
void Call_neighbor(){
int threadsPerBlock =128 ;
int blocksPerGrid =(N+threadsPerBlock-1)/threadsPerBlock ;
//size_t size;
neighbor<<<blocksPerGrid, threadsPerBlock>>>(d_U,d_UR,d_UL);

}
__global__ void neighbor(float *U,float *UR,float *UL){

int i= blockDim.x * blockIdx.x +threadIdx.x ;

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
