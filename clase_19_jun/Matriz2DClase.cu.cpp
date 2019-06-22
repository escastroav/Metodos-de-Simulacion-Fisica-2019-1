// Programa que hace c=a+b para a,b,c vectores en CUDA
#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Lx 128
#define Ly 128
#define Nx 8
#define Ny 8
#define Q 5
const int Mx=(Lx+Nx-1)/Nx;
const int My=(Ly+Ny-1)/Nx;

const double W0=1.0/3;

const double C=0.5; // C<0.707 celdas/click
const double TresC2=3*C*C;
const double AUX0=1-TresC2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//--------------------KERNELS----------------
__constant__ float d_w[Q];
__constant__ int d_Vx[Q];
__constant__ int d_Vy[Q];
__global__ void IncrementarMatriz(float *d_a,size_t pitcha){
  int ix,iy; float *a;
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;

  a=d_a+(ix*pitcha)/sizeof(float)+iy;
  
  (*a)++;
}
//------------------- CLASES ----------------
class LatticeBoltzmann{
private:
  float h_w[Q]; int h_Vx[Q], h_Wy[Q];
  
  float h_f0[Lx][Ly]; float*d_f0; size_t pitchf0; 
  float h_f1[Lx][Ly]; float*d_f1; size_t pitchf1; 
  float h_f2[Lx][Ly]; float*d_f2; size_t pitchf2; 
  float h_f3[Lx][Ly]; float*d_f3; size_t pitchf3; 
  float h_f4[Lx][Ly]; float*d_f4; size_t pitchf4;

  float h_f0new[Lx][Ly]; float*d_f0new; size_t pitchf0new; 
  float h_f1new[Lx][Ly]; float*d_f1new; size_t pitchf1new; 
  float h_f2new[Lx][Ly]; float*d_f2new; size_t pitchf2new; 
  float h_f3new[Lx][Ly]; float*d_f3new; size_t pitchf3new; 
  float h_f4new[Lx][Ly]; float*d_f4new; size_t pitchf4new;
 
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  void Inicie(void);
  void Incremente(void);
  void Imprimase(void);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  h_w[0]=W0; h_w[1]=h_w[2]=h_w[3]=h_w[4]=(1-W0)/4.0;
  cudaMemcpyToSymbol(d_w,h_w,Q*sizeof(float),0,cudaMemcpyHostToDevice); //Enviarlos Al Device
  //Cargar los vectores
  h_Vx[0]=0;
  h_Vy[0]=0;
  
  h_Vx[1]=1;  h_Vx[2]=0;  h_Vx[3]=-1; h_Vx[4]=0;
  h_Vy[1]=0;  h_Vy[2]=1;  h_Vy[3]=0;  h_Vy[4]=-1;
  
  cudaMemcpyToSymbol(d_Vx,h_Vx,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  
  //Construir las matrices virtuales en el Device
  
  cudaMallocPitch((void**) &d_f0,&pitchf0,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f1,&pitchf1,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f2,&pitchf2,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f3,&pitchf3,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f4,&pitchf4,Ly*sizeof(float),Lx);

  cudaMallocPitch((void**) &d_f0new,&pitchf0new,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f1new,&pitchf1new,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f2new,&pitchf2new,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f3new,&pitchf3new,Ly*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f4new,&pitchf4new,Ly*sizeof(float),Lx);
  
}
LatticeBoltzmann::~LatticeBoltzmann(void){
  cudaFree(d_f0);
  cudaFree(d_f1);
  cudaFree(d_f2);
  cudaFree(d_f3);
  cudaFree(d_f4);
  
  cudaFree(d_f0new);
  cudaFree(d_f1new);
  cudaFree(d_f2new);
  cudaFree(d_f3new);
  cudaFree(d_f4new);
}
void LatticeBoltzmann::Inicie(void){
  //Llevar al Device
  cudaMemcpy2D(d_a,pitcha,h_a,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
}
void LatticeBoltzmann::Imprimase(void){
  int ix,iy;
  //Devolver los datos al Host
  cudaMemcpy2D(h_a,Ly*sizeof(float),d_a,pitcha,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  //Mostrar
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++)
      cout<<h_a[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;
}
void LatticeBoltzmann::Incremente(void){
  //Procesar en el Device
  dim3 ThreadsPerBlock(Nx,Ny,1);
  dim3 BlocksPerGrid(Mx,My,1);
  IncrementarMatriz<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,pitcha);
}
// ----------------- FUNCIONES GLOBALES ----------
int main(void){
  LatticeBoltzmann Ondas;

  Ondas.Inicie();
  Ondas.Imprimase();
  Ondas.Incremente();
  Ondas.Imprimase();

  return 0;
}
