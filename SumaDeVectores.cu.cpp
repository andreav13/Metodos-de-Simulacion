#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>

using namespace std;

#define Lx 128  //Numero total de hilos
#define Nx 64   //Numero de hilos por bloque
const int Mx=(Lx+Nx-1)/Nx;   //Numero de bloques

__global__ void SumeDosVectores(float *d_a, float *d_b, float *d_c){
  int ix;
  ix=blockIdx.x*blockDim.x+threadIdx.x;
}


int main(void){

  //DECLARAR LAS MATRICES
  float h_a[Lx];   float *d_a; cudaMalloc((void**) &d_a, Lx*sizeof(float));
  float h_b[Lx];   float *d_b; cudaMalloc((void**) &d_b, Lx*sizeof(float));
  float h_c[Lx];   float *d_c; cudaMalloc((void**) &d_c, Lx*sizeof(float));
  int ix;

  //INICIALIZAR LOS DATOS
  //Cargar en el host
  for(ix=0;ix<Lx;ix++) {
    h_a[ix]=ix;
    h_b[ix]=2*ix;
  }
  //Enviarlo al device
  //        a donde, de donde,   tamano,        flag    
  cudaMemcpy(d_a, h_a, Lx*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, h_b, Lx*sizeof(float), cudaMemcpyHostToDevice);
    
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  SumeDosVectores<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,d_b,d_c);

  //IMPRIMIR LOS DATOS
  //Devolverlos al host
  cudaMemcpy(h_test,d_test,Lx*sizeof(float),cudaMemcpyDeviceToHost);
  //Imprimirlos
  for(ix=0;ix<Lx;ix++) cout<<h_test[ix]<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_test);

  
  return 0;

}
