#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>

using namespace std;

#define Lx 16
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx; //Define la cantidad de bloques suficientes para cubrir la grilla.

//--------------KERNELS----------------
__global__ void Sume(float *d_a, float *d_b, float *d_c)//Debe ser siempre de tipo void
{
  int ix;

  ix = blockIdx.x*blockDim.x + threadIdx.x; //Recorriendo en un arreglo unidimensional agrupados por bloque de 8.

  d_c[ix] = d_a[ix] + d_b[ix];
}

int main(void)
{
  //DECLARAR LAS MATRICES
  float h_a[Lx],h_b[Lx],h_c[Lx];
  float *d_a; cudaMalloc((void**) &d_a,Lx*sizeof(float));
  float *d_b; cudaMalloc((void**) &d_b,Lx*sizeof(float));
  float *d_c; cudaMalloc((void**) &d_c,Lx*sizeof(float));

  //INICIALIZAR DATOS
  //cargar en el host
  for(int ix=0; ix<Lx;ix++)
    {
      h_a[ix] = ix; h_b[ix] = 2*ix; 
    }
  //Enviarlos al device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_b,h_b,Lx*sizeof(float),cudaMemcpyHostToDevice);
  
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlockPerGrid(Mx,1,1);
  Sume<<<BlockPerGrid,ThreadsPerBlock>>>(d_a,d_b,d_c);

  //DEVOLVER AL HOST
  cudaMemcpy(h_c,d_c,Lx*sizeof(float),cudaMemcpyDeviceToHost);

  //Imprimir
  for(int ix=0; ix < Lx; ix++)
    cout << h_c[ix] << endl;

  //Liberar Memoria
  cudaFree(d_a);
  cudaFree(d_b);
  cudaFree(d_c);
  
  return 0;
}
