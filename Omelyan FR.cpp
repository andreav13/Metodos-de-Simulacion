// Métodos de Integración.
#include<iostream>
#include<fstream>
#include<cmath>
#include "Vector.h"
using namespace std;

const double G=1.0;
const int N=2;
const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;
class Cuerpo;
class Colisionador;

//---------------------Clase Cuerpo--------------------------_//

class Cuerpo{                   // punto y coma obligatorio despues del corchete
private:
  vector3D r,V,F;   double m,R;
public:
  void  Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0);
  void BorrarFuerza(void);
  void AgregueFuerza(vector3D F0);
  void Mueva_r(double dt,double Constante);
  void Mueva_v(double dt,double Constante);
  void Dibujese(void);
  double Getx(void){return r.x();};        // Función Inline
  double Gety(void){return r.y();};
  friend class Colisionador;
};

void  Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0){
  r.cargue(x0,y0,z0); V.cargue(Vx0,Vy0,Vz0); m=m0; R=R0;
}

void Cuerpo::BorrarFuerza(void){
  F.cargue(0,0,0);
}

void Cuerpo::AgregueFuerza(vector3D F0){
  F+=F0;
}

void Cuerpo::Mueva_r(double dt,double Constante){
  r+=V*(Constante*dt);
}

void Cuerpo::Mueva_v(double dt,double Constante){
  V+=F*(Constante*dt)/m;
}

void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}


//---------------------Clase Colisionador--------------------------_//

class Colisionador{
private:
public:
void  CalculeTodasLasFuerzas(Cuerpo* Planeta);
void  CalculeLaFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2);
};


// El * es para permitir pasar a la función un array de planetas.
void Colisionador:: CalculeTodasLasFuerzas(Cuerpo* Planeta){
  int i,j;
  // Borrar todas las fuerzas
  for(i=0;i<N;i++){Planeta[i].BorrarFuerza();}
  // Calculas todas las fuerzas entre parejas de planetas.
  for(i=0;i<N;i++){
    for(j=i+1;j<N;j++){
      CalculeLaFuerzaEntre(Planeta[i],Planeta[j]);
    }
  }
}


void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2){
  vector3D F1,dr=Planeta2.r-Planeta1.r;
   double aux=G*Planeta1.m*Planeta2.m*pow(norma2(dr),-1.5);
   F1=dr*aux;
   Planeta1.AgregueFuerza(F1); Planeta2.AgregueFuerza(F1*(-1));
}

//---------------------Funciones Globales--------------------------_//

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'MiPlanetaVector.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-120:120]"<<endl;
  cout<<"set yrange [-120:120]"<<endl;
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
  double t,dt=50.0;
  int Ndibujos,tdibujo;
  Cuerpo Planeta[N]; int i;
  Colisionador Newton;
  
  double m0=10, m1=1, r=1000;
  double R0=10,R1=5;
  
  double M=m0+m1;
  double x0=r*(-m1/M), x1=x0+r;
  
  double omega=sqrt(G*M/(r*r*r)),Vy0=omega*x0,Vy1=omega*x1,T=2*M_PI/omega,tmax=1.1*T;
  
  //  InicieAnimacion();
  // Ndibujos=500;
  //--------------------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0);
  Planeta[0].Inicie(x0,0,0,0,Vy0,0,m0,R0);
  Planeta[1].Inicie(x1,0,0,0,Vy1,0,m1,R1);

  for (t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    
    // if (tdibujo>tmax/Ndibujos){
    //   InicieCuadro();
    //   Planeta.Dibujese();
    //   TermineCuadro();
    //   tdibujo=0;
    // }
    
    cout<<Planeta[1].Getx()<<" "<<Planeta[1].Gety()<<endl;
    
    //Muevase con Omelyan FR.

    for(i=0;i<N;i++){ Planeta[i].Mueva_r(dt,Zeta);}
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_v(dt,(1-2*Lambda)/2);}
    for(i=0;i<N;i++){Planeta[i].Mueva_r(dt,Xi);}
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_v(dt,Lambda);}
    for(i=0;i<N;i++){Planeta[i].Mueva_r(dt,1-2*(Xi+Zeta));}
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_v(dt,Lambda);}
    for(i=0;i<N;i++){Planeta[i].Mueva_r(dt,Xi);}
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_v(dt,(1-2*Lambda)/2);}
    for(i=0;i<N;i++){Planeta[i].Mueva_r(dt,Zeta);}
  }
  
  return 0;
}



// Algoritmo predictor corrector
