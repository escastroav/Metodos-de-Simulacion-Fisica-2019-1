#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>

using namespace std;

#define Lx 16
#define Ly 8
#define Nx 8
#define Ny 8
const int Mx=(Lx+Nx-1)/Nx;
const int My=(Ly+Ny-1)/Ny; //Define la cantidad de bloques suficientes para cubrir la grilla.

//--------------KERNELS----------------
__global__ void SumeleUno(float *d_f, size_t pitchf)//Debe ser siempre de tipo void
{
  //Identificar el procesador que hara el calculo.
  int ix,iy;
  ix = blockIdx.x*blockDim.x + threadIdx.x; //Recorriendo en un arreglo unidimensional agrupados por bloque de 8.
  iy = blockIdx.y*blockDim.y + threadIdx.y;
  //Cuál es la direccion de memoria del elemento
  float *aux;
  aux=d_f+(ix*pitchf)/sizeof(float)+iy; //aux es &(d_f[ix][iy])
  //Usarla
  (*aux)++; 
}

int main(void)
{
  //DECLARAR LAS MATRICES
  float h_f[Lx][Ly];
  float *d_f; size_t pitchf;//Tamaño del paquete de BUS para enviar la información de forma organizada.
  cudaMallocPitch((void**) &d_f,&pitchf,Lx*sizeof(float),Lx);
  
  //INICIALIZAR DATOS
  //cargar en el host
  for(int ix=0; ix<Lx;ix++)
    for(int iy=0; iy<Ly;iy++)
      h_f[ix][iy]=ix*Ly+iy;

   //Imprimir
  for(int ix=0; ix < Lx; ix++)
    {
      for(int iy=0; iy < Ly; iy++)
	cout << h_f[ix][iy] << " ";
      cout << endl;
    }cout << endl;
  
  //Enviarlos al device
  cudaMemcpy2D(d_f, pitchf, h_f, Ly*sizeof(float), Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,Ny,1);
  dim3 BlockPerGrid(Mx,My,1);
  SumeleUno<<<BlockPerGrid,ThreadsPerBlock>>>(d_f,pitchf);

  //DEVOLVER AL HOST
  cudaMemcpy2D(h_f,Ly*sizeof(float),d_f,pitchf,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);

  //Imprimir
  for(int ix=0; ix < Lx; ix++)
    {
      for(int iy=0; iy < Ly; iy++)
	cout << h_f[ix][iy] << " ";
      cout << endl;
    }cout << endl;
  
  //Liberar Memoria
  cudaFree(d_f);
  
  return 0;
}
