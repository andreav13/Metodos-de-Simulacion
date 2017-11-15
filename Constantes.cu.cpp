// Difusion 1D 
#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

//--------------------KERNELS----------------
__constant__ float d_w[5];
__global__ void Test(float *d_test){
  (*d_test)=d_w[0];
}

int main(void){
  float h_w[5]; // w[i]
  /*Test*/ float h_test[1]; float*d_test; cudaMalloc((void**) &d_test,sizeof(float));

  //Cargar las constantes en el Host
  h_w[0]=1.0/3;    h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6;
  //Enviarlas al Device
  cudaMemcpyToSymbol(d_w,h_w,5*sizeof(float),0,cudaMemcpyHostToDevice);

  //TEST
  //Llevar al Device
  /*Test*/ h_test[0]=0; cudaMemcpy(d_test,h_test,sizeof(float),cudaMemcpyHostToDevice);
  
  //Cambiar en el Device
  dim3 ThreadsPerBlock(1,1,1);
  dim3 BlocksPerGrid(1,1,1);
  Test<<<BlocksPerGrid,ThreadsPerBlock>>>(d_test);

  //Devolver al Host
  cudaMemcpy(h_test,d_test,sizeof(float),cudaMemcpyDeviceToHost);

  cout<<h_test[0]<<endl;

  cudaFree(d_test);

  return 0;
}
