#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<sys/time.h>


#define High 100.0
#define Long 100.0
#define wide 100.0
#define NX 2560
#define NY 1
#define NZ 1
#define N  (NX*NY*NZ)
#define DX (High/NX)
#define DY  (Long/NY)
#define DZ  (wide/NZ)
#define DT  (0.05*DX)
#define NO_STEPS 2000
#define g   9.81
//#define GAMA 1.4
//#define R 287.0
//#define Cv (R/(GAMA-1))
//#define Cp (Cv+R)


#ifndef EXTERNAL
#define EXTERNAL extern
#endif


void Send_To_Device();
void Get_From_Device();
void Free_Memory();
void Allocate_Memory();
void Call_Unew();
void Call_neighbor();

EXTERNAL float *h_P;
EXTERNAL float *h_U;
EXTERNAL float *d_U;
EXTERNAL float *d_UR;
EXTERNAL float *d_UL;



