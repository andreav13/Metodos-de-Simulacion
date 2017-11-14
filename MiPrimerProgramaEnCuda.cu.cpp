#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>

using namespace std;

__constant__ float d_w[5];
__global__ void SumeleUno(float *d_test){
  (*d_test)+=d_w[0];
}


int main(void){

  //DECLARAR LAS MATRICES
  float h_w[5]; //w[i]
  
  float h_test[1];    //h por host, d por device
  float *d_test; cudaMalloc((void**) &d_test, sizeof(float));

  //INICIALIZAR LOS DATOS
  //Cargar en el host
  h_test[0]=10;
  h_w[0]=1./3; h_w[1]=h_w[2]=h_w[3]=h_w[4]=1./6;
  //Enviarlo al device
  //        a donde, de donde,   tamano,        flag    
  cudaMemcpy(d_test, h_test, sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_w,h_w,5*sizeof(float),0,cudaMemcpyHostToDevice);
    
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(1,1,1);
  dim3 BlocksPerGrid(1,1,1);
  SumeleUno<<<BlocksPerGrid,ThreadsPerBlock>>>(d_test); //funcionAQUIEN(parametro)

  //IMPRIMIR LOS DATOS
  //Devolverlos al host
  cudaMemcpy(h_test,d_test,sizeof(float),cudaMemcpyDeviceToHost);
  //Imprimirlos
  cout<<h_test[0]<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_test);

  
  return 0;

}
