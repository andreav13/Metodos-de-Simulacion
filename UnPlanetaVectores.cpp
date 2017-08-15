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
  void Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0);
  void CalculeFuerza(void);
  void Muevase(double dt);
  void Dibujese(void);
  double Getx(void){return x;};
  double Gety(void){return y;};
};


//Funciones de la clase cuerpo
void Cuerpo::Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0){
  x=x0; y=y0; Vx=Vx0; Vy=Vy0; m=m0; R=R0;
}


void Cuerpo::CalculeFuerza(void){
  double r2=x*x+y*y;
  double aux=GM*m*pow(r2,-1.5);
  Fx=-aux*x; Fy=-aux*y;
}


void Cuerpo::Muevase(double dt){
  x+=Vx*dt; y+=Vy*dt;
  Vx+=Fx*dt/m; Vy+=Fy*dt/m;
}


void Cuerpo::Dibujese(void){
  cout<<", "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
}


//Funciones globales
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'MiPlaneta.gif'"<<endl;
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

  double t, tdibujo, tmax, dt=0.05;
  int Ndibujos=1000;
  Cuerpo Planeta;

  double m=1, R=5;
  double r, omega, V, T;
  r=100; omega=sqrt(GM*pow(r,-3)); V=omega*r; T=2*M_PI/omega; tmax=1.1*T;


  //InicieAnimacion();
  
  //(x0, y0, Vx0, Vy0, m0, R0)
  Planeta.Inicie(r, 0, 0, 0.5*V, m, R);
  
  for(t=0, tdibujo=0;t<tmax;t+=dt, tdibujo+=dt){
    /*if(tdibujo>tmax/Ndibujos){
      InicieCuadro();
      Planeta.Dibujese();
      TermineCuadro();
      tdibujo=0;
      }*/

    //cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
    Planeta.CalculeFuerza();
    Planeta.Muevase(dt);
  }
  
  return 0;

}


//PARA CORRER:    ./a.out | gnuplot
