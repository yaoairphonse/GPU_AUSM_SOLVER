#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<sys/time.h>


#define High 1.0
#define Long 1.0
#define wide 1.0
#define NX 100
#define NY 100
#define NZ 100
#define N  (NX*NY*NZ)
#define DX (High/NX)
#define DY  DX
#define DZ  DY
#define DT 5e-07
#define NO_STEPS 8000
#define GAMA 1.4
#define R 287.0
#define Cv (R/(GAMA-1))
#define Cp (Cv+R)


#ifndef EXTERNAL
#define EXTERNAL extern
#endif


void Send_To_Device();
void Get_From_Device();
void Free_Memory();
void Allocate_Memory();
void Call_Unew();
void Call_neighbor();

EXTERNAL float *h_body;
EXTERNAL float *h_P;
EXTERNAL float *h_U;


//EXTERNAL float *d_P;
EXTERNAL float *d_U;

EXTERNAL float *h_UR ;
EXTERNAL float *h_UL ;
EXTERNAL float *h_UT ;
EXTERNAL float *h_UD ;
EXTERNAL float *h_UF ;
EXTERNAL float *h_UB ;

EXTERNAL  float *d_body;
EXTERNAL  float *d_UR ;
EXTERNAL  float *d_UL ;
EXTERNAL  float *d_UT ;//top
EXTERNAL  float *d_UD ;//down
EXTERNAL  float *d_UF ;
EXTERNAL  float *d_UB ;
