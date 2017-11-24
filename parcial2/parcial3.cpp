#include <iostream>
#include <fstream>
#include <cmath> 
#include "Random64.h"

using namespace std;

const int Lx=2000;
const int Ly=1;
const int Q=5;
const double W0=1/3.;
const double C=0.5;  
const double tresC2=3*C*C;
const double AUX0=1-tresC2*(1-W0);
const double tau=0.5;
const double Utau=1./tau;
const double UmUtau=1-Utau;


class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q]; //V[alpha][i] alpha=1 es x, alpha=0 es y
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; //f[ix][iy][i]
  double rho_min[Lx], rho_max[Lx];
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, bool UseNew, double sigma);
  double Jx(int ix, int iy);
  double Jy(int ix, int iy);
  double fequilibrio(int i, double rho0, double Jx0, double Jy0);
  void Inicie(double rho0, double Jx0, double Jy0);
  double GetSigma(int ix, int iy, int t);
  void Colisione(int t);
  void Adveccione(void);
  void ActualiceEnvolventes(int t);
  void ImprimaRho(void);
  void ImprimaA(void);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  w[0]=W0;
  w[1]=w[2]=w[3]=w[4]=1/6.;

  V[0][0]=0;
  V[1][0]=0;

  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1;  V[0][4]=0;
  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;   V[1][4]=-1;

  for (int ix=0; ix<Lx; ix++){
    rho_min[ix]=0; rho_max[ix]=0;
  }
  
}

double LatticeBoltzmann::rho(int ix, int iy, bool UseNew, double sigma){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
    suma+=fnew[ix][iy][i];
    else
      suma+=f[ix][iy][i];
  return suma+1/.2*sigma;
}

double LatticeBoltzmann::Jx(int ix, int iy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    suma+=V[0][i]*f[ix][iy][i];
  return suma;
}

double LatticeBoltzmann::Jy(int ix, int iy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    suma+=V[1][i]*f[ix][iy][i];
  return suma;
}

double LatticeBoltzmann::fequilibrio(int i, double rho0, double Jx0, double Jy0){
  if(i==0)
    return AUX0*rho0;
  else
    return w[i]*(tresC2*rho0+3*(V[0][i]*Jx0+V[1][i]*Jy0));
  
}

void LatticeBoltzmann::Inicie(double rho0, double Jx0, double Jy0){
int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){
	f[ix][iy][i]=fequilibrio(i,rho0,Jx0,Jy0);
      }
}

double LatticeBoltzmann::GetSigma(int ix, int iy, int t){
  double A=1, lambda=1000, omega=2*M_PI*C/lambda;
  if(ix==0)
    return A*sin(omega*t);
  else
    return 0;
}

void LatticeBoltzmann::Colisione(int t){ //de f a fnew
  int ix,iy,i,i_intercambio=0; double rho0,Jx0,Jy0,sigma,D=0.6;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){ //para cada celda

      if(ix==Lx-1)
	for(i=0;i<Q;i++){
	  if(i==1 or i==2) i_intercambio=i+2; if(i==3 or i==4) i_intercambio=i-2;
	  fnew[ix][iy][i]=D*f[ix][iy][i_intercambio]; 
	}

      else{
	sigma=GetSigma(ix,iy,t);
	rho0=rho(ix,iy,false,sigma); Jx0=Jx(ix,iy); Jy0=Jy(ix,iy);
	for(i=0;i<Q;i++) //para cada direccion
	  fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*fequilibrio(i,rho0,Jx0,Jy0);
      }
    }
}
   
void LatticeBoltzmann::Adveccione(void){ //de fnew a f
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}

void LatticeBoltzmann::ActualiceEnvolventes(int t){
  for(int ix=1;ix<Lx-1;ix++){
    int iy=0;
    double sigma = GetSigma(ix,iy,t); double rho0=rho(ix,iy,true,sigma);
    if (rho0<rho_min[ix]) rho_min[ix]=rho0;
    if (rho0>rho_max[ix]) rho_max[ix]=rho0;
  }
}

void LatticeBoltzmann::ImprimaRho(void){
  ofstream RhoMax("rhomax.dat");
  ofstream RhoMin("rhomin.dat");
  for(int ix=1;ix<Lx-1;ix++){
    RhoMax<<ix<<" "<<rho_max[ix]<<endl;
    RhoMin<<ix<<" "<<rho_min[ix]<<endl;
  }
  RhoMax.close();
  RhoMin.close();
}

void LatticeBoltzmann::ImprimaA(void){
  double Amax=0, Amin=30, SWR, Cabs; 
  for(int ix=1;ix<Lx-1;ix++){
    if (rho_max[ix]>Amax) Amax=rho_max[ix];
    if (rho_max[ix]<Amin) Amin=rho_max[ix];
  }

  SWR=Amax/Amin;
  Cabs=4*SWR/((SWR+1)*(SWR+1));
  
  cout<<"SWR es "<< SWR << " y Cabsorcion es "<< Cabs << endl;
}


//-----------------------Funciones Globales-----------------------


int main(void){

  LatticeBoltzmann Ondas;
  int t,tmax=80000;

  double rho0=0,Jx0=0,Jy0=0;
  
  //Inicie
  Ondas.Inicie(rho0,Jx0,Jy0);

  //Corra
  for(t=0;t<tmax;t++){
    Ondas.Colisione(t);
    Ondas.Adveccione();
  }
  
  for(t=tmax;t<tmax+2200;t++){
    Ondas.Colisione(t);
    Ondas.Adveccione();
    Ondas.ActualiceEnvolventes(t);
  }
  
  Ondas.ImprimaRho();
  Ondas.ImprimaA();

  
  return 0;
}
