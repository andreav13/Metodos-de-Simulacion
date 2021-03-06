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

__constant__ float d_w[5];
__device__ float SumeleUno(float x){
  return x+1;
}
__global__ void SumeleAlgo(float *d_test){
  (*d_test)+=SumeleUno(d_w[0]);
}


int main(void){

  //DECLARAR LAS MATRICES
  float h_w[5]; //w[i]
  
  float h_test[Lx];    //h por host, d por device
  float *d_test; cudaMalloc((void**) &d_test, sizeof(float));
  int ix;

  //INICIALIZAR LOS DATOS
  //Cargar en el host
  for(ix=0;ix<Lx;ix++) h_test[ix]=ix;
  h_test[0]=10;
  h_w[0]=1./3; h_w[1]=h_w[2]=h_w[3]=h_w[4]=1./6;
  //Enviarlo al device
  //        a donde, de donde,   tamano,        flag    
  cudaMemcpy(d_test, h_test, Lx*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_w,h_w,5*sizeof(float),0,cudaMemcpyHostToDevice);
    
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  SumeleAlgo<<<BlocksPerGrid,ThreadsPerBlock>>>(d_test); //funcionAQUIEN(parametro)

  //IMPRIMIR LOS DATOS
  //Devolverlos al host
  cudaMemcpy(h_test,d_test,Lx*sizeof(float),cudaMemcpyDeviceToHost);
  //Imprimirlos
  for(ix=0;ix<Lx;ix++) cout<<h_test[ix]<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_test);

  
  return 0;

}
