#include <iostream>
#include <fstream> //save to file
#include <cmath> //math functions
#include "Vector.h"

using namespace std;

const double GM=1;
class Cuerpo;

//Clase Cuerpo
class Cuerpo{
private:
  vector3D r,V,F;
  double m,R;
  
public:
  void Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0);
  void CalculeFuerza(void);
  void Muevase(double dt);
  void Dibujese(void);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
};


//Funciones de la clase cuerpo
void Cuerpo::Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0){
  r.cargue(x0,y0,z0);
  V.cargue(Vx0,Vy0,Vz0);
  m=m0; R=R0;
}


void Cuerpo::CalculeFuerza(void){
  double aux=GM*m*pow(norma2(r),-1.5);
  F=(-aux)*r;
}


void Cuerpo::Muevase(double dt){
  r+=V*dt;
  V+=F*(dt/m);
}


void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}


//Funciones globales
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl;
  //cout<<"set output 'MiPlaneta.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-120:120]"<<endl;
  cout<<"set yrange[-120:120]"<<endl;
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

  double t, tdibujo, tmax, dt=0.01;
  //int Ndibujos=1000;
  Cuerpo Planeta;

  double m=1, R=5;
  double r, omega, V, T;
  r=100; omega=sqrt(GM*pow(r,-3)); V=omega*r; T=2*M_PI/omega; tmax=1.1*T;


  //InicieAnimacion();
  
  //(x0, y0, z0, Vx0, Vy0, Vz0, m0, R0)
  Planeta.Inicie(r, 0, 0, 0, 0.5*V, 0, m, R);
  
  for(t=0, tdibujo=0;t<tmax;t+=dt, tdibujo+=dt){
    /*if(tdibujo>tmax/Ndibujos){
      InicieCuadro();
      Planeta.Dibujese();
      TermineCuadro();
      tdibujo=0;
      }*/

    cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
    Planeta.CalculeFuerza();
    Planeta.Muevase(dt);
  }
  
  return 0;

}
