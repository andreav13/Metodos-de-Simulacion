#include <iostream>
#include <fstream> //save to file
#include <cmath> //math functions
#include "Vector.h"
#include "Random64.h"

using namespace std;

const double Lx=100, Ly=100;
const int Nx=1, Ny=1, N=Nx*Ny;
const double E=1.0, r0=10;

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;
  
class Cuerpo;
class Colisionador;

//Clase Cuerpo
class Cuerpo{
private:
  vector3D r,V,F;
  double m,R,re;
  
public:
  void Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0);
  void BorreFuerza(void);
  void AgregueFuerza(double F0);
  void Mueva_r(double dt, double Constante);
  void Mueva_V(double dt, double Constante);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  double Getr(void){return re;};
  void Dibujese(void);

  friend class Colisionador;
};


//Clase Colisionador
class Colisionador{
private:

public:
  void CalculeTodasLasFuerzas(Cuerpo* Grano);
};
  
//Funciones de la clase cuerpo
void Cuerpo::Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0){
  r.cargue(x0,y0,z0);
  re=norma(r);
  V.cargue(Vx0,Vy0,Vz0);
  m=m0; R=R0;
}


void Cuerpo::BorreFuerza(void){
  F.cargue(0,0,0);
}


void Cuerpo::AgregueFuerza(double F0){
  vector3D r_unitario=r/re;
  F+=F0*r_unitario;
}


void Cuerpo::Mueva_r(double dt, double Constante){
  r+=V*(Constante*dt);
}


void Cuerpo::Mueva_V(double dt, double Constante){
  V+=F*(Constante*dt)/m;
}


void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}


void Colisionador::CalculeTodasLasFuerzas(Cuerpo* Grano){
  int i,j;  
  
  for(i=0;i<N;i++){
    Grano[i].BorreFuerza();
  }
  //Agregue la fuerza de Lennard-Jones
  for(i=0;i<N;i++){
    Grano[i].AgregueFuerza(12*E/Grano[i].re*(pow((r0/Grano[i].re),12)-pow((r0/Grano[i].re),6)));
  }
  
}


int main(void){
  
  int i,j;
  double t, dt=0.1;
  Cuerpo Grano[N];
  Colisionador Newton;

  double kbT=0.5;
  double m0=1.0, R0=2.5, x0=10, y0=0, V0=pow(2*kbT/m0,0.5), Vx0=V0, Vy0=0;

  double tmax=100;

  double KbT;
  for(KbT=0.025;KbT<0.5;KbT+=0.05){
    
  //   Inicializa granos
  //           (x0, y0, z0, Vx0, Vy0, Vz0, m0, R0)
  for(i=0;i<Nx;i++)
    for(j=0;j<Ny;j++){
      Grano[i+Nx*j].Inicie(x0, y0, 0, pow(2*KbT/m0,0.5), Vy0, 0, m0, R0);
      //Grano[i+Nx*j].Inicie(x0, y0, 0, Vx0, Vy0, 0, m0, R0);
    }

  double rmin=0, rmax=0, rcentral;
  for (t=0;t<tmax;t+=dt){

    //cout<<t<<" "<<Grano[0].Getx()<<endl;

    if(Grano[0].Getx()>rmax) rmax=Grano[0].Getx();
    if(Grano[0].Getx()<rmin) rmin=Grano[0].Getx();

    rcentral=(rmax+rmin)/2.;

      
  //Muevase con Omelyan FR.
    for(i=0;i<N;i++){
      Grano[i].Mueva_r(dt, Zeta);
    }
    Newton.CalculeTodasLasFuerzas(Grano);
    
    for(i=0;i<N;i++){
      Grano[i].Mueva_V(dt, (1-2*Lambda)/2);
    }
    for(i=0;i<N;i++){
      Grano[i].Mueva_r(dt, Xi);
    }
   
    Newton.CalculeTodasLasFuerzas(Grano);
    
    for(i=0;i<N;i++){
      Grano[i].Mueva_V(dt, Lambda);
    }
    for(i=0;i<N;i++){
      Grano[i].Mueva_r(dt, 1-2*(Xi+Zeta)); 
    }
    Newton.CalculeTodasLasFuerzas(Grano);
    
    for(i=0;i<N;i++){
      Grano[i].Mueva_V(dt, Lambda);
    }
    for(i=0;i<N;i++){
      Grano[i].Mueva_r(dt, Xi);
    }
    
    Newton.CalculeTodasLasFuerzas(Grano);
    
    for(i=0;i<N;i++){
      Grano[i].Mueva_V(dt, (1-2*Lambda)/2);
    }
    for(i=0;i<N;i++){
      Grano[i].Mueva_r(dt, Zeta);
    }

    

  }
  cout<<KbT<<" "<<rcentral<<endl;

  }
  
  return 0;

}
