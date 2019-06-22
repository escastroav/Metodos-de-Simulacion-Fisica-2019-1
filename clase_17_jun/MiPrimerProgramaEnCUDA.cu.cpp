#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>

using namespace std;
//--------------KERNEL----------------
__global__ void Test(float *d_test)//Debe ser siempre de tipo void
{
  (*d_test)++;
}

int main(void)
{
  //DECLARAR LAS MATRICES
  float h_test[1];
  float *d_test; cudaMalloc((void**) &d_test,sizeof(float));

  //INICIALIZAR DATOS
  //cargar en el host
  h_test[0]=3.4;
  //Enviarlos al device
  cudaMemcpy(d_test,h_test,sizeof(float),cudaMemcpyHostToDevice);//cudaMemcpyHostToDevice es una constante que puede ser 0 o 1

  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(1,1,1);
  dim3 BlockPerGrid(1,1,1);
  Test<<<BlockPerGrid,ThreadsPerBlock>>>(d_test);

  //DEVOLVER AL HOST
  cudaMemcpy(h_test,d_test,sizeof(float),cudaMemcpyDeviceToHost);

  //Imprimir
  cout << h_test[0] << endl;

  //Liberar Memoria
  cudaFree(d_test);
  
  return 0;
}
