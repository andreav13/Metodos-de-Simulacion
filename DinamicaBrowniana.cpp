#include <iostream>
#include <fstream> 
#include <cmath> 
#include "Random64.h"

using namespace std;

const double Deltat=0.01, L=100, T=300, MASA= 22.8916, e=1, D=1.132, kb=0.826;//unidades: picosegundos, Amstromg, Kelvin, u.m.a.
const double Gamma=kb*T/(MASA*D), sigma=sqrt(2*D*Deltat), dtU2mGamma=Deltat/(2*MASA*Gamma);

class Cuerpo;

//Clase Cuerpo
class Cuerpo{
private:
  double x,Fx,Fxold,m,R;
  
public:
  void Inicie(double x0, double m0, double R0);
  void CalculeFuerza(double E);
  void Dibujese(void);
  void Muevase(Crandom & ran64);
  double Getx(void){return x;}; 
};


//Funciones de la clase cuerpo
void Cuerpo::Inicie(double x0, double m0, double R0){
  x=x0; m=m0; R=R0; Fx=0;
}


void Cuerpo::CalculeFuerza(double E){
  Fxold=Fx; Fx=e*E;
}


void Cuerpo::Muevase(Crandom & ran64){
  x+=dtU2mGamma*(3*Fx-Fxold) + ran64.gauss(0,sigma); 
}

void Cuerpo::Dibujese(void){
  cout<<", "<<x<<"+"<<R<<"*cos(t),"<<0<<"+"<<R<<"*sin(t)";
}


//Funciones globales
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'DinBrowniana.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:110]"<<endl;
  cout<<"set yrange[-10:10]"<<endl;
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


//Main
int main(void){

  double tdibujo;  int Ndibujos=1000;
  double t,  tmax=40;
  Crandom ran64(1);
  double E=0;
  Cuerpo Na;//Sodio

  InicieAnimacion();
  
  //(x0, m0, R0)
  Na.Inicie(L/2, MASA, 4);
  Na.CalculeFuerza(E);

  for(t=0, tdibujo=0;t<tmax;t+=Deltat, tdibujo+=Deltat){
    if(tdibujo>tmax/Ndibujos){
      InicieCuadro();
      Na.Dibujese();
      TermineCuadro();
      tdibujo=0;

      Na.CalculeFuerza(E);
      Na.Muevase(ran64);
      }

    
  }
  
  return 0;

}
