#include <iostream>
#include <fstream> //save to file
#include <cmath> //math functions
#include "Vector.h"
#include "Random64.h"

using namespace std;

const double K=1e4, Gamma=50;
const double g=9.8;
const double Lx=100, Ly=100;
const int Nx=1, Ny=1, N=Nx*Ny;
  

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;
  
class Cuerpo;
class Colisionador;

//Clase Cuerpo
class Cuerpo{
private:
  vector3D r,V,F,omega,tau;
  double m,R,theta,I;
  
public:
  void Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0,double theta0, double omega0, double m0, double R0);
  void BorreFuerzayTorque(void);
  void AgregueFuerza(vector3D F0);
  void AgregueTorque(vector3D tau0);
  void Mueva_r(double dt, double Constante);
  void Mueva_V(double dt, double Constante);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  void Dibujese(void);

  friend class Colisionador;
};


//Clase Colisionador
class Colisionador{
private:

public:
  void CalculeTodasLasFuerzas(Cuerpo* Grano);
  void CalculeLaFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2);
};
  
//Funciones de la clase cuerpo
void Cuerpo::Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0,double theta0, double omega0, double m0, double R0){
  r.cargue(x0,y0,z0);
  V.cargue(Vx0,Vy0,Vz0);
  omega.cargue(0,0,omega0);
  theta=theta0;
  m=m0; R=R0; I=2/5.*m*R*R;
}


void Cuerpo::BorreFuerzayTorque(void){
  F.cargue(0,0,0);
  tau.cargue(0,0,0);
}


void Cuerpo::AgregueFuerza(vector3D F0){
  F+=F0;
}

void Cuerpo::AgregueTorque(vector3D tau0){
  tau+=tau0;
}

void Cuerpo::Mueva_r(double dt, double Constante){
  r+=V*(Constante*dt);
  theta+=omega.z()*Constante*dt;
}


void Cuerpo::Mueva_V(double dt, double Constante){
  V+=F*(Constante*dt)/m;
  omega+=tau*Constante*dt/I;
}


void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}


void Colisionador::CalculeTodasLasFuerzas(Cuerpo* Grano){
  int i,j;
  vector3D g_vector; g_vector.cargue(0,-g,0);
  for(i=0;i<N;i++){
    Grano[i].BorreFuerzayTorque();
  }
  //Agregue la fuerza de la gravedad
  for(i=0;i<N;i++){
    Grano[i].AgregueFuerza((Grano[i].m)*g_vector);
  }
  
  //Calcular todas las fuerzas entre parejas de planetas
  for(i=0;i<N;i++){
    for(j=i+1;j<N+4;j++){
      CalculeLaFuerzaEntre(Grano[i], Grano[j]);
    }
  }
}


void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2){
  vector3D F2,Fn,runitario,n,Vc,Vcn,Vct,t,r21 = Grano2.r-Grano1.r;
  double s,m1,m2,R1,R2,m12,componenteVcn,normaVct,componenteFn;
  double ERFF=1e-8,d21=norma(r21);
  s=Grano1.R+Grano2.R-d21;

  if(s>0){

    //Geometria y dinamica del contacto
    m1=Grano1.m; m2=Grano2.m; m12=(m1*m2)/(m1+m2);
    R1=Grano1.R; R2=Grano2.R;
    n=r21/d21;

    //Calcular velocidad de contacto y el vector tangente
    Vc=(Grano2.V-Grano1.V); //velocidad del punto de contacto
                            // -> dividimos en componentes normal y tangencial
    componenteVcn=Vc*n;Vcn=n*componenteVcn; Vct=Vc-Vcn; normaVct=norma(Vct);
    if(normaVct<ERFF) t.cargue(0,0,0);
    else t=Vct/normaVct;

    //FUERZAS NORMALES
    //Fuerza de Hertz
    componenteFn=K*pow(s,1.5);
    //Disipación Plástica
    componenteFn-=m12*sqrt(s)*Gamma*componenteVcn;
    if(componenteFn<0) componenteFn=0;
    Fn=n*componenteFn;

    
    
    F2=Fn;
    Grano1.AgregueFuerza(F2*(-1)); Grano2.AgregueFuerza(F2);
  }
}


void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'granos_que_giran.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-10:110]"<<endl;
  cout<<"set yrange [-10:110]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;
}

void InicieCuadro(void){
  cout<<"plot 0,0 ";
  cout<<" , "<<100.1/7<<"*t,0";
  cout<<" , "<<100.1/7<<"*t,100";
  cout<<" , 0,"<<100.1/7<<"*t";
  cout<<" , 100,"<<100.1/7<<"*t";
}
void TermineCuadro(void){
  cout<<endl;
}



int main(void){
  
  int i,j;
  double t, dt=1e-3;
  int Ndibujos;
  double tdibujo;
  Cuerpo Grano[N+4];
  Colisionador Newton;
  Crandom ran64(1);

  double m0=1, R0=3, V=10, omega0=10;
  double Rpared=10000, Mpared=1000;

  double T=Lx/V, tmax=5*T;

  double dx=Lx/(Nx+1), dy=Ly/(Ny+1), theta;

  InicieAnimacion();
  Ndibujos=500;

  //   Paredes
  //           (x0, y0, z0, Vx0, Vy0, Vz0, theta0, omega0, m0, R0)

  //pared arriba
  Grano[N].Inicie(Lx/2, Ly+Rpared, 0, 0, 0, 0, 0, 0, Mpared, Rpared);
  //pared abajo
  Grano[N+1].Inicie(Lx/2, -Rpared, 0, 0, 0, 0, 0, 0, Mpared, Rpared);
  //pared derecha
  Grano[N+2].Inicie(Lx+Rpared, Ly/2, 0, 0, 0, 0, 0, 0, Mpared, Rpared);
  //pared izquierda
  Grano[N+3].Inicie(-Rpared, Ly/2, 0, 0, 0, 0, 0, 0, Mpared, Rpared);

  
  //   Inicializa granos
  //           (x0, y0, z0, Vx0, Vy0, Vz0, theta0, omega0, m0, R0)
  for(i=0;i<Nx;i++)
    for(j=0;j<Ny;j++){
      theta=2*M_PI*ran64.r();
      Grano[i+Nx*j].Inicie((i+1)*dx, (j+1)*dy, 0, 0, 0, 0, 0, omega0, m0, R0);
    }

  
  for (t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    
    if (tdibujo>tmax/Ndibujos){
      InicieCuadro();
      for(i=0;i<N;i++){Grano[i].Dibujese();}
      TermineCuadro();
      tdibujo=0;
    }
    
  
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
  
  return 0;

}
