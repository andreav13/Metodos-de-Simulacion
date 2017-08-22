#include <iostream>
#include <fstream> //save to file
#include <cmath> //math functions
#include "Vector.h"

using namespace std;

const double GM=1;
const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;
  
class Cuerpo;

//Clase Cuerpo
class Cuerpo{
private:
  vector3D r,V,F;
  double m,R;
  
public:
  void Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0);
  void CalculeFuerza(void);
  void Mueva_r(double dt, double Constante);
  void Mueva_V(double dt, double Constante);
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


void Cuerpo::Mueva_r(double dt, double Constante){
  r+=V*(Constante*dt);
}


void Cuerpo::Mueva_V(double dt, double Constante){
  V+=F*(Constante*dt);
}



int main(void){

  double t, tdibujo, tmax, dt=15;
  Cuerpo Planeta;

  double m=1, R=5;
  double r, omega, V, T;
  r=100; omega=sqrt(GM*pow(r,-3)); V=omega*r; T=2*M_PI/omega; tmax=1.1*T;


  
  //         (x0, y0, z0, Vx0, Vy0, Vz0, m0, R0)
  Planeta.Inicie(r, 0, 0, 0, 0.5*V, 0, m, R);
  
  for(t=0, tdibujo=0;t<tmax;t+=dt, tdibujo+=dt){
    cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
    
    Planeta.Mueva_r(dt, Zeta);
    Planeta.CalculeFuerza();
    Planeta.Mueva_V(dt, (1-2*Lambda)/2);
    Planeta.Mueva_r(dt, Xi);
    Planeta.CalculeFuerza();
    Planeta.Mueva_V(dt, Lambda);
    Planeta.Mueva_r(dt, 1-2*(Xi+Zeta));
    Planeta.CalculeFuerza();
    Planeta.Mueva_V(dt, Lambda);
    Planeta.Mueva_r(dt, Xi);
    Planeta.CalculeFuerza();
    Planeta.Mueva_V(dt, (1-2*Lambda)/2);
    Planeta.Mueva_r(dt, Zeta);
  }
  
  return 0;

}
