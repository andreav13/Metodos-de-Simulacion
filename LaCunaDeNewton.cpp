#include<iostream>
#include<fstream>
#include<cmath>
#include "Vector.h"
using namespace std;

const double g=9.8;
const int N=6;
const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;
class Cuerpo;
class Colisionador;

//---------------------Clase Cuerpo--------------------------_//

class Cuerpo{
private:
  //vector3D tau, omega;
  double theta, omega, tau, m, R, L, I, xcorrido;
public:
  void  Inicie(double theta0, double omega0, double m0, double R0, double L0, double x0corrido);
  void InicieFuerza(void);
  void BorrarFuerza(void);
  void AgregueFuerza(double tau0);
  void Mueva_theta(double dt,double Constante);
  void Mueva_omega(double dt,double Constante);
  void Dibujese(void);
  //double Getx(void){return r.x();};    
  //double GetFx(void);
  friend class Colisionador;
};

void  Cuerpo::Inicie(double theta0, double omega0, double m0, double R0, double L0, double x0corrido){
  omega=omega0;
  theta=theta0;
  m=m0;
  R=R0;
  L=L0;
  I=m*L*L;
  xcorrido=x0corrido;
}

void Cuerpo::InicieFuerza(){
  tau=-m*g*L*sin(theta);
}

void Cuerpo::BorrarFuerza(void){
  tau=0;
}

void Cuerpo::AgregueFuerza(double tau0){
  tau+=tau0;
}

void Cuerpo::Mueva_theta(double dt,double Constante){
  theta+=omega*(Constante*dt);
}

void Cuerpo::Mueva_omega(double dt,double Constante){
  omega+=tau*(Constante*dt)/I;
}

void Cuerpo::Dibujese(void){
  cout<<", "<<xcorrido+L*sin(theta)<<"+"<<R<<"*cos(t),"<<-L*cos(theta)<<"+"<<R<<"*sin(t) , "
      <<xcorrido<<"+"<<(xcorrido+L*sin(theta)-xcorrido)/7.0<<"*t,"<<0<<"-"<<L*cos(theta)/7.0<<"*t";
  
}


//---------------------Clase Colisionador--------------------------_//

class Colisionador{
private:
public:
  void  CalculeTodasLasFuerzas(Cuerpo* Pendulo);
  void  CalculeLaFuerzaEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2);
};


// El * es para permitir pasar a la funciÃ³n un array de planetas.
void Colisionador:: CalculeTodasLasFuerzas(Cuerpo* Pendulo){
  int i,j;
  // Borrar todas las fuerzas
  for(i=0;i<N;i++){Pendulo[i].InicieFuerza();}
  //Agregar fuerzas externas
  //for(i=0;i<N;i++){Pendulo[i].AgregueFuerza(-Pendulo[i].m*g*sin(Pendulo[i].theta));}
  // Calculas todas las fuerzas entre parejas de péndulos.
  /*for(i=0;i<N;i++){
    for(j=i+1;j<N;j++){
      CalculeLaFuerzaEntre(Pendulo[i],Pendulo[j]);
    }
    }*/
}


void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Pendulo1, Cuerpo & Pendulo2){
  /*vector3D F1,dr=Pendulo2.r-Pendulo1.r;
  double aux=G*Pendulo1.m*Pendulo2.m*pow(norma2(dr),-1.5);
  F1=dr*aux;
  Pendulo1.AgregueFuerza(F1); Pendulo2.AgregueFuerza(F1*(-1));*/
}

//---------------------Funciones Globales--------------------------_//

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'MiPenduloVector.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-12:22]"<<endl;
  cout<<"set yrange [-12:0]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;
}

void InicieCuadro(void){
  cout<<"plot 0,0 ";
}
void TermineCuadro(void){
  cout<<endl;
}

//-----------------Programa Principal----------------//

int main(void){
  double t,dt=1.0;
  int Ndibujos,tdibujo;
  Cuerpo Pendulo[N]; int i;
  Colisionador Newton;
  
  double m0=1, L0=10,R0=1, x0corrido=0, theta0=-15*M_PI/180;
  double T=2*M_PI*sqrt(L0/g),tmax=3*T;
  
  InicieAnimacion();
  Ndibujos=500;
  
  //-------------(theta0,omega0,m0,R0,L0,x0corrido);
  Pendulo[0].Inicie(theta0,0,m0,R0,L0,x0corrido);
  for(i=1;i<N;i++)Pendulo[i].Inicie(0,0,m0,R0,L0,2*R0*i);
  

  for (t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    
    if (tdibujo>tmax/Ndibujos){
      InicieCuadro();
      for(i=0;i<N;i++){Pendulo[i].Dibujese();}
      TermineCuadro();
      tdibujo=0;
    }
    
    
    //Muevase con Omelyan FR.

    for(i=0;i<N;i++){ Pendulo[i].Mueva_theta(dt,Zeta);}
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0;i<N;i++){Pendulo[i].Mueva_omega(dt,(1-2*Lambda)/2);}
    for(i=0;i<N;i++){Pendulo[i].Mueva_theta(dt,Xi);}
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0;i<N;i++){Pendulo[i].Mueva_omega(dt,Lambda);}
    for(i=0;i<N;i++){Pendulo[i].Mueva_theta(dt,1-2*(Xi+Zeta));}
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0;i<N;i++){Pendulo[i].Mueva_omega(dt,Lambda);}
    for(i=0;i<N;i++){Pendulo[i].Mueva_theta(dt,Xi);}
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0;i<N;i++){Pendulo[i].Mueva_omega(dt,(1-2*Lambda)/2);}
    for(i=0;i<N;i++){Pendulo[i].Mueva_theta(dt,Zeta);}
  }
  
  return 0;
}



// Algoritmo predictor corrector
