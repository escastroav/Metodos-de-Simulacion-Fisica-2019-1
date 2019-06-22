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
const int Mx=(Lx+Nx-1)/Nx;
const int My=(Ly+Ny-1)/Nx;

#define Q 5
#define C (0.5)
#define TresC2 (0.75)
#define AUX0 (0.5)

#define tau (0.5)
#define Utau (2.0)
#define UmUtau (-1.0)


//--------------------KERNELS----------------
__constant__ float d_w[Q];
__constant__ int d_Vx[Q];
__constant__ int d_Vy[Q];

__device__ float d_rho(float f0,float f1,float f2,float f3,float f4){
  return f0+f1+f2+f3+f4;
}

__device__ float J_x(float f0,float f1,float f2,float f3,float f4){
  return f0*d_Vx[0]+f1*d_Vx[1]+f2*d_Vx[2]+f3*d_Vx[3]+f4*d_Vx[4];
}

__device__ float J_y(float f0,float f1,float f2,float f3,float f4){
  return f0*d_Vy[0]+f1*d_Vy[1]+f2*d_Vy[2]+f3*d_Vy[3]+f4*d_Vy[4];
}

__device__ float d_feq0(float rho0, float Jx0,float Jy0,int i){
  return rho0*AUX0;
}
__device__ float d_feqi(float rho0, float Jx0,float Jy0,int i){
    return d_w[i]*(TresC2*rho0+3*(d_Vx[i]*Jx0+d_Vy[i]*Jy0));
}

//kernels de evolucion
__global__ void d_Colisione(float *f0, size_t pitchf0,
			  float *f1, size_t pitchf1,
			  float *f2, size_t pitchf2,
			  float *f3, size_t pitchf3,
			  float *f4, size_t pitchf4,      
			  float *f0new,size_t pitchf0new,			  
			  float *f1new,size_t pitchf1new,			  
			  float *f2new,size_t pitchf2new,			  
			  float *f3new,size_t pitchf3new,			  
			  float *f4new,size_t pitchf4new,
			  float rhoFuente)
{
  //Quien soy yo?
  int ix, iy; //float *f0,*f1,*f2,*f3,*f4,*f0new,*f1new,*f2new,*f3new,*f4new;
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;
  f0=d_f0+(ix*pitchf0)/sizeof(float)+iy; f0new=d_fnew0+(ix*pitchf0new)/sizeof(float)+iy;
  f1=d_f1+(ix*pitchf1)/sizeof(float)+iy; f1new=d_fnew1+(ix*pitchf1new)/sizeof(float)+iy;
  f2=d_f2+(ix*pitchf2)/sizeof(float)+iy; f2new=d_fnew2+(ix*pitchf2new)/sizeof(float)+iy;
  f3=d_f3+(ix*pitchf3)/sizeof(float)+iy; f3new=d_fnew3+(ix*pitchf3new)/sizeof(float)+iy;
  f4=d_f4+(ix*pitchf4)/sizeof(float)+iy; f4new=d_fnew4+(ix*pitchf4new)/sizeof(float)+iy;
  //Procesar Datos
  //calcular cantidades macroscopicas
  float rho0,Jx0,Jy0;
  rho0=d_rho(*f0,*f1,*f2,*f3,*f4); Jx0=d_Jx(*f0,*f1,*f2,*f3,*f4); Jy0=d_Jy(*f0,*f1,*f2,*f3,*f4);
  //evolucionar
  (*f0new)=UmUtau*(*f0)+Utau*d_feq0(rho0,Jx0,Jy0,0);
  (*f1new)=UmUtau*(*f1)+Utau*d_feqi(rho0,Jx0,Jy0,1);
  (*f2new)=UmUtau*(*f2)+Utau*d_feqi(rho0,Jx0,Jy0,2);
  (*f3new)=UmUtau*(*f3)+Utau*d_feqi(rho0,Jx0,Jy0,3);
 (*f4new)=UmUtau*(*f4)+Utau*d_feqi(rho0,Jx0,Jy0,4); 
}
__global__ void d_ImponerCampos(float *f0, size_t pitchf0,
				float *f1, size_t pitchf1,
				float *f2, size_t pitchf2,
				float *f3, size_t pitchf3,
				float *f4, size_t pitchf4,      
				float *f0new,size_t pitchf0new,			  
				float *f1new,size_t pitchf1new,			  
				float *f2new,size_t pitchf2new,			  
				float *f3new,size_t pitchf3new,			  
				float *f4new,size_t pitchf4new,
				float rhoFuente)
{
  //Quien soy yo?
  int ix, iy; //float *f0,*f1,*f2,*f3,*f4,*f0new,*f1new,*f2new,*f3new,*f4new;
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;
  //Procesar Datos
  float rho0,Jx0,Jy0;
  if(ix==Lx/2 && iy==Ly/2)
    {
      f0=d_f0+(ix*pitchf0)/sizeof(float)+iy; f0new=d_fnew0+(ix*pitchf0new)/sizeof(float)+iy;
      f1=d_f1+(ix*pitchf1)/sizeof(float)+iy; f1new=d_fnew1+(ix*pitchf1new)/sizeof(float)+iy;
      f2=d_f2+(ix*pitchf2)/sizeof(float)+iy; f2new=d_fnew2+(ix*pitchf2new)/sizeof(float)+iy;
      f3=d_f3+(ix*pitchf3)/sizeof(float)+iy; f3new=d_fnew3+(ix*pitchf3new)/sizeof(float)+iy;
      f4=d_f4+(ix*pitchf4)/sizeof(float)+iy; f4new=d_fnew4+(ix*pitchf4new)/sizeof(float)+iy;
      //Procesar Datos
      
      //Calcular cantidades macroscópicas
      rho0=rhoFuente; Jx0=d_Jx(*f0,*f1,*f2,*f3,*f4); Jy0=d_Jy(*f0,*f1,*f2,*f3,*f4);
    }
  //calcular cantidades macroscopicas
  
}  
__global__ void IncrementarMatriz(float *d_a,size_t pitcha){
  int ix,iy; float *a;
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;

  a=d_a+(ix*pitcha)/sizeof(float)+iy;
  
  (*a)++;
}
//------------------- CLASES ----------------
class LatticeBoltzmann{
private:
  float h_w[Q]; int h_Vx[Q], h_Vy[Q];
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
  double h_rho(int ix,int iy,bool UseNew);
  double h_Jx(int ix,int iy,bool UseNew);
  double h_Jx(int ix,int iy,bool UseNew);
  void h_Colisione(void);
  double h_feq(double rho0,double Jx0,double Jy0,int i);
  void Inicie(void);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //construir los pesos
  //Cargar los pesos en el host
  h_w[0]=1.0/3; h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6;
  cudaMemcpyToSymbol(d_w,h_w,Q*sizeof(float),0,cudaMemcpyHostToDevice); //enviamos al device
  
  //Cargar los vectores
  h_Vx[0]=0;
  h_Vy[0]=0;

  h_Vx[1]=1;  h_Vx[2]=0;  h_Vx[3]=-1; h_Vx[4]=0;
  h_Vy[1]=0;  h_Vy[2]=1;  h_Vy[3]=0;  h_Vy[4]=-1;
  cudaMemcpyToSymbol(d_Vx,h_Vx,Q*sizeof(int),0,cudaMemcpyHostToDevice); //enviarlos al device
  cudaMemcpyToSymbol(d_Vy,h_Vy,Q*sizeof(int),0,cudaMemcpyHostToDevice); //enviarlos al device
  
  //Construir las matrices f en el Device
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
  cudaFree(d_f0); cudaFree(d_f1); cudaFree(d_f2); cudaFree(d_f3); cudaFree(d_f4);
cudaFree(d_f0new); cudaFree(d_f1new); cudaFree(d_f2new); cudaFree(d_f3new); cudaFree(d_f4new);
}

double LatticeBoltzmann::h_rho(int ix,int iy,bool UseNew){
    if(UseNew)
      return h_f0new[ix][iy][i]+h_f1new[ix][iy][i]+h_f2new[ix][iy][i]+h_f3new[ix][iy][i]+h_f4new[ix][iy][i];
    else
      return h_f0[ix][iy][i]+h_f1[ix][iy][i]+h_f2[ix][iy][i]+h_f3[ix][iy][i]+h_f4[ix][iy][i];
}

double LatticeBoltzmann::h_Jx(int ix,int iy,bool UseNew){
  if(UseNew)
    return h_Vx[0]*h_f0new[ix][iy]+h_Vx[1]*h_f1new[ix][iy]+h_Vx[2]*h_f2new[ix][iy]+h_Vx[3]*h_f3new[ix][iy]+h_Vx[4]*h_f4new[ix][iy];

  else
    return h_f0[ix][iy]*h_Vx[0]+h_f1[ix][iy]*h_Vx[1]+h_f2[ix][iy]*h_Vx[2]+h_f3[ix][iy]*h_Vx[3]+h_f4[ix][iy]*h_Vx[4];
}
double LatticeBoltzmann::h_Jy(int ix,int iy,bool UseNew){
  if(UseNew)
    return h_Vy[0]*h_f0new[ix][iy]+h_Vy[1]*h_f1new[ix][iy]+h_Vy[2]*h_f2new[ix][iy]+h_Vy[3]*h_f3new[ix][iy]+h_Vy[4]*h_f4new[ix][iy];

  else
    return h_f0[ix][iy]*h_Vy[0]+h_f1[ix][iy]*h_Vy[1]+h_f2[ix][iy]*h_Vy[2]+h_f3[ix][iy]*h_Vy[3]+h_f4[ix][iy]*h_Vy[4];
}
double LatticeBoltzmann::h_feq(double rho0,double Jx0,double Jy0,int i){
  if(i==0)
    return rho0*AUX0;
  else
    return h_w[i]*(TresC2*rho0+3*(h_Vx[i]*Jx0+h_Vy[i]*Jy0));
}




void LatticeBoltzmann::Inicie(void){
    for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      h_f0[ix][iy]=feq(rho0,Jx0,Jy0,0);
      h_f1[ix][iy]=feq(rho0,Jx0,Jy0,1);
      h_f2[ix][iy]=feq(rho0,Jx0,Jy0,2);
      h_f3[ix][iy]=feq(rho0,Jx0,Jy0,3);
      h_f4[ix][iy]=feq(rho0,Jx0,Jy0,4);
    }
    cudaMemcpy2D(d_f0,pitchf0,h_f0,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
    cudaMemcpy2D(d_f1,pitchf1,h_f1,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
    cudaMemcpy2D(d_f2,pitchf2,h_f2,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
    cudaMemcpy2D(d_f3,pitchf3,h_f3,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
    cudaMemcpy2D(d_f4,pitchf0,h_f4,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);

}

void LatticeBoltzmann::h_Colisione(void){
  //Procesar en el Drive
  dim3 ThreadsPerBlock(Nx,Ny,1);
  dim3 BlocksPerGrid(Mx,My,1);
  d_Colisione<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f0,pitchf0,
						 d_f1,pitchf1,
						 d_f2,pitchf2,
						 d_f3,pitchf3,
						 d_f4,pitchf4,
						 d_f0new,pitchf0new,
						 d_f1new,pitchf1new,
						 d_f2new,pitchf2new,
						 d_f3new,pitchf3new,
						 d_f4new,pitchf4new)
}
void LatticeBoltzmann::ImponerCampos(int t){
  double A, lambda,omega;
  A=10; lambda=10; omega=2*M_PI/lambda;
  d_ImponerCampos<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f0,pitchf0,
						 d_f1,pitchf1,
						 d_f2,pitchf2,
						 d_f3,pitchf3,
						 d_f4,pitchf4,
						 d_f0new,pitchf0new,
						 d_f1new,pitchf1new,
						 d_f2new,pitchf2new,
						 d_f3new,pitchf3new,
						 d_f4new,pitchf4new,
						     A*sin(t*omega))
  
}

// ----------------- FUNCIONES GLOBALES ----------
int main(void){
  LatticeBoltzmann Ondas;
  /*
  Ondas.Inicie();
  Ondas.Imprimase();
  Ondas.Incremente();
  Ondas.Imprimase();
  */
  return 0;
}
