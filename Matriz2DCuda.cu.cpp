// Programa que hace c=a+b para a,b,c vectores en CUDA
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
const int My=(Ly+Ny-1)/Ny;

//--------------------KERNELS----------------
__device__ float SumeleUno(float x){
  return x+1;
}

__global__ void SumarConstante(float *d_a, size_t pitcha){
  int ix,iy; float *aux;
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;
  
  aux=d_a+(ix*pitcha)/sizeof(float)+iy;

  (*aux)=SumeleUno((*aux));  // (*aux) es &d_a[ix][iy]
}

int main(void){
  int ix,iy;
  //DECLARAR
  //Declarar matrices en el Host
  float h_a[Lx][Ly]; 
  //Declarar matrices en el Device
  float*d_a; size_t pitcha; cudaMallocPitch((void**) &d_a,&pitcha,Ly*sizeof(float),Lx);

  //INICIALIZAR
  //Cargar valores en el Host
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      h_a[ix][iy]=Ly*ix+iy;

  //Mostrar
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++)
      cout<<h_a[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;
  
  //Llevar al Device
  cudaMemcpy2D(d_a,pitcha,h_a,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  
  //PROCESAR
  //Procesar en el Device
  dim3 ThreadsPerBlock(Nx,Ny,1);
  dim3 BlocksPerGrid(Mx,My,1);
  SumarConstante<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,pitcha);

  //MOSTRAR RESULTADOS
  //Devolver al Host
  cudaMemcpy2D(h_a,Ly*sizeof(float),d_a,ptcha,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      cout<<h_a[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;

  cudaFree(d_a);  

  return 0;
}
