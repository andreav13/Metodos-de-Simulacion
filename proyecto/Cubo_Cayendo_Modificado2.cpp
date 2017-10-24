#include <iostream>
#include <fstream>
#include <cmath> 
#include "Vector.h"
#include "Random64.h"

using namespace std;

const double A=10., L=100., R=1.;
const double Gamma=0.50, Kcundall=10,K=1e4;
const double g=9.8;
const int N=27;

// ------------Elastic Moduli.  1/E = (1-v1²)/E1 + (1-v2²)/E2-------//

//const double v1=0.27,v2=0.27,E1=200,E2=200;         // Objeto y superficie de Metal.
//const double v1=0.33,v2=0.33,E1=117,E2=117;         // Objeto y superficie de Cobre.
const double v1=0.4999,v2=0.4999,E1=0.01,E2=0.01;     // Objeto y superficie de Caucho.

const double UnoSobreE=(1-v1*v1)/E1 + (1-v2*v2)/E2;
const double E=1/UnoSobreE;

//-----------Constantes Para corres Omelyan----------------//

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;

//-------------Definición de Clases----------------------//

class Cuerpo;
class Colisionador;

//Clase Cuerpo
class Cuerpo{
private:
  vector3D r,rot,a,b,c,V,F,omega,tau;
  double m,theta,I;
  
public:
  void Inicie(double x0, double y0, double z0, double phi0, double alpha0, double beta0, double a10, double a20, double a30, double b10, double b20, double b30,  double c10, double c20, double c30, double Vx0, double Vy0, double Vz0,double theta0, double omega0, double m0, double I0);
  void BorreFuerzayTorque(void);
  void AgregueFuerza(vector3D F0);
  void AgregueTorque(vector3D tau0);
  void Mueva_r(double dt, double Constante);
  void Mueva_V(double dt, double Constante);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  double Getz(void){return r.z();};
  double Getaz(void){return a.z();};
  double Getbz(void){return b.z();};
  double Getcz(void){return c.z();};
  double GetVx(void){return V.x();};
  double GetVy(void){return V.y();};
  double GetVz(void){return V.z();};
  double GetFz(void){return F.z();};
  void Dibujar_esfera(void);
  void Dibujar_cilindro(void);
  void Dibujar_plano(void);
  friend class Colisionador;
};


//Clase Colisionador
class Colisionador{
private:
public:
  void Inicie(void);
  void CalculeTodasLasFuerzas(Cuerpo* Grano, double dt);
  void CalculeLaFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2);
  void CalculeLaFuerzaEntrePlanos(Cuerpo * Grano, double dt);
  };
  
//Funciones de la clase cuerpo
void Cuerpo::Inicie(double x0, double y0, double z0, double phi0, double alpha0, double beta0, double a10, double a20, double a30, double b10, double b20, double b30,  double c10, double c20, double c30, double Vx0, double Vy0, double Vz0,double theta0, double omega0, double m0, double I0){
  r.cargue(x0,y0,z0);
  rot.cargue(phi0,alpha0,beta0);
  a.cargue(a10,a20,a30);
  b.cargue(b10,b20,b30);
  c.cargue(c10,c20,c30);
  V.cargue(Vx0,Vy0,Vz0);
  theta=theta0;
  omega.cargue(0,0,omega0);
  m=m0; I=I0;
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
  a+=V*(Constante*dt);
  b+=V*(Constante*dt);
  c+=V*(Constante*dt);
  theta+=omega.z()*Constante*dt;
}


void Cuerpo::Mueva_V(double dt, double Constante){
  V+=F*(Constante*dt)/m;
  //  omega+=tau*Constante*dt/I;
}


void Cuerpo::Dibujar_esfera(void){
  cout<<", "<<r.x()<<"+"<<R<<"*sin(u)*cos(v),"<<r.y()<<"+"<<R<<"*sin(u)*sin(v),"<<r.z()<<"+"<<R<<"*cos(u)";
}

void Cuerpo::Dibujar_cilindro(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(v)*cos("<<rot.x()<<")*cos("<<rot.y()<<")"<<"+t(u)*sin("<<rot.x()<<")*cos("<<rot.y()<<")-"<<R<<"*sin(v)*sin("<<rot.y()<<"),"
      <<r.y()<<"+"<<R<<"*cos(v)*cos("<<rot.x()<<")*sin("<<rot.y()<<")*cos("<<rot.z()<<")"<<"+t(u)*sin("<<rot.x()<<")*sin("<<rot.y()<<")*cos("<<rot.z()<<")+"<<R<<"*sin(v)*cos("<<rot.y()<<")*cos("<<rot.z()<<")+"<<R<<"*cos(v)*sin("<<rot.x()<<")*sin("<<rot.z()<<")-cos("<<rot.x()<<")*t(u)*sin("<<rot.z()<<"),"
      <<r.z()<<"+"<<R<<"*cos(v)*cos("<<rot.x()<<")*sin("<<rot.y()<<")*sin("<<rot.z()<<")+sin("<<rot.x()<<")*t(u)*sin("<<rot.y()<<")*sin("<<rot.z()<<")+"<<R<<"*sin(v)*cos("<<rot.y()<<")*sin("<<rot.z()<<")-"<<R<<"*cos(v)*sin("<<rot.x()<<")*cos("<<rot.z()<<")+cos("<<rot.x()<<")*t(u)*cos("<<rot.z()<<")";
}

void Cuerpo::Dibujar_plano(void){
  cout<<", "<<a.x()<<"+m(u)*("<<b.x()<<"-"<<a.x()<<")+n(v)*("<<c.x()<<"-"<<a.x()<<"),"
      <<a.y()<<"+m(u)*("<<b.y()<<"-"<<a.y()<<")+n(v)*("<<c.y()<<"-"<<a.y()<<"),"
      <<a.z()<<"+m(u)*("<<b.z()<<"-"<<a.z()<<")+n(v)*("<<c.z()<<"-"<<a.z()<<")";
}


void Colisionador::CalculeTodasLasFuerzas(Cuerpo* Grano,double dt){
  int i,j;
  vector3D g_vector; g_vector.cargue(0,0,g);
  
  for(i=0;i<N;i++){Grano[i].BorreFuerzayTorque();}
  //Agregue la fuerza de la gravedad.
  for(i=0;i<N-1;i++){Grano[i].AgregueFuerza(-(Grano[i].m)*g_vector);}
  //Agregue la fuerza viscosa.
  for(i=0;i<N;i++){Grano[i].AgregueFuerza(-Gamma*Grano[i].V);}
  // Calcule fuerza entre el objeto y el plano.
  CalculeLaFuerzaEntrePlanos(Grano,dt);
  
}

void Colisionador::CalculeLaFuerzaEntrePlanos(Cuerpo * Grano,double dt){
  vector3D F2,Fn,n;
  double componenteFn;

  //Vector unitario normal a la superficie.
  n.cargue(0,0,-1);   

  //-------------------------------------------------------------------------------
  //Colisión entre las caras del cubo y la superficie. Tipo de Colisión Plano-Plano
  //-------------------------------------------------------------------------------
  
  for (int i = 20;i<N;i++){
    if (Grano[i].Getaz()<-L/2 and Grano[i].Getbz()<-L/2 and Grano[i].Getcz()<-L/2){
      double s_plano=-Grano[i].Getaz()-L/2;
      
      //componenteFn=E*L*s_plano;                        //Ecuación para la Colisión plano-plano?? dada por el profe.
                                                         // No es lo suficientemente fuerte y la fuerza de gravedad+viscosa predominan
                                                         // permitiendo que la tierra se trague al cubito.
      componenteFn=K*pow(s_plano,1.5);
      if(componenteFn<0){ componenteFn=0;}
      Fn=n*componenteFn;
      F2=Fn;
    }
  }

  //--------------------------------------------------------------------------------------------------------
  //Colisión entre las esferas (vertices) del cubo y la superficie(esfera). Tipo de Colisión Esfera-Esfera.;
  //--------------------------------------------------------------------------------------------------------
  
  // cout<<Grano[6].Getx()<<" "<<Grano[6].Gety()<<" "<<Grano[6].Getz()<<endl;
  
  for (int i = 0;i<8;i++){
    double r1,r2=-L/2;
    if(Grano[i].Getx()<-L/2 or Grano[i].Gety()<-L/2 or Grano[i].Getz()<-L/2){        //Condición de contacto Esfera-Superficie.
      if(Grano[i].Getx()<-L/2){r1=Grano[i].Getx();}                                  //Determinar que parte de la esfera hizo contacto.
      else if(Grano[i].Gety()<-L/2){r1=Grano[i].Gety();}
      else if(Grano[i].Getz()<-L/2){r1=Grano[i].Getz();}
      double Rpared=10000;                                                          //Necesaria para colisión esfera-esfera.
      double d_esfera= -r1+r2;                                                      //Distancia que penetra la esfera en el plano(otra esfera)
      double s_esfera=(Rpared + R)-d_esfera;                                        //Es raro pero se requiere para la colisión segun el profe
      double R_efectivo=(Rpared*R)/(Rpared+R);                                      // Necesario para el calculo de K
      double K_esfera=(4/3)*E*pow(R_efectivo,0.5);                                  // K calculado manualmente
      componenteFn=K_esfera*pow(s_esfera,1.5);
      if(componenteFn<0){ componenteFn=0;}
      Fn=n*componenteFn;
      F2=Fn;
    }
  }

  //---------------------------------------------------------------------------------------------------------------------------
  // Colisión entre los cilindros (lados) del cubo y la superficie (cilindro).Tipo de Colisión Cilindro-Cilindro-Ejes Paralelos.
  //---------------------------------------------------------------------------------------------------------------------------
  
  //for (int i = 8;i<20;i++){
  
  //Condición de contacto cilindro-Superficie, supongo que debe ser cuando la coordenada z de la base del cilindro este por debajo de -L/2
  //Fn=(M_PI/4)*E*L*s       con   L=longitud del cilindro.
  
    // }
  
  for(int i=0;i<N-1;i++){Grano[i].AgregueFuerza(F2*(-1));}
}


void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'cubo_modificado_con_Fuerzas.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-60:60]"<<endl;
  cout<<"set yrange [-60:60]"<<endl;
  cout<<"set zrange [-60:12]"<<endl;
  cout<<"set size 1.1,1.1"<<endl;                           // El tamaño de la imagen aumenta, es como hacer un poco de zoom.
  cout<<"set size ratio 1"<<endl;                           // 1 o -1, No se observan cambios importantes
  cout<<"set view equal xyz"<<endl;
  cout<<"set hidden3d front"<<endl;
  
  //cout<<"set pm3d depthorder interpolate 0,0"<<endl;      // Esta en el Tetraedro, se puede comentar! No se observan cambios importantes.
  //cout<<"set palette rgb 3,3,3"<<endl;                    // Esta en el Tetraedro, se puede comentar! No se observan cambios importantes.
  
  cout<<"set parametric"<<endl;
  cout<<"set urange [0:pi]"<<endl;
  cout<<"set vrange [0:2*pi]"<<endl;
  cout<<"t(u)=10*u/pi"<<endl;
  cout<<"m(u)=u/pi"<<endl;
  cout<<"n(v)=v/(2*pi)"<<endl;
  cout<<"set isosamples 8,8"<<endl;                         //20 o 8,8, Al cambiar a 20 el cuerpo se colorea como
                                                            //solidos pero compila mas lento, se ve mejor con 8,8. 
}

void InicieCuadro(void){
  cout<<"splot 0,0,0 ";
}

void TermineCuadro(void){
  cout<<endl;
}

int main(void){

  double t, dt=1e-2, tmax=20;
  int Ndibujos;
  double tdibujo;
  Cuerpo parte[N];
  Colisionador Newton;
  Crandom ran64(1);

  double me=1,mc=1,mp=1,mp_grande=1000;
  double Ie=2/5.*me*R*R, Ic=1/2.*mc*R*R, Ip=1/6.*mp*A*A, Ip_grande=1/6.*mp_grande*L*L;
  double V=0, Vy=0, theta0=M_PI, omega0=10;

  Ndibujos=100;
  InicieAnimacion();
  

  //Inicio de cuerpos
  //x0, y0, z0,   phi0, alpha0, beta0,
  //a10,   a20,   a30,   b10,   b20,   b30,    c10,   c20,   c30,
  //Vx0,   Vy0,   Vz0,  theta0,   omega0,   m0,  I

  //ESFERAS
  parte[0].Inicie(0,0,0, 0,0,0, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,me,Ie);
  parte[1].Inicie(A,0,0, 0,0,0, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,me,Ie);
  parte[2].Inicie(0,A,0, 0,0,0, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,me,Ie);
  parte[3].Inicie(0,0,A, 0,0,0, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,me,Ie);
  parte[4].Inicie(0,A,A, 0,0,0, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,me,Ie);
  parte[5].Inicie(A,0,A, 0,0,0, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,me,Ie);
  parte[6].Inicie(A,A,0, 0,0,0, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,me,Ie);
  parte[7].Inicie(A,A,A, 0,0,0, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,me,Ie);
  
  //CILINDROS
  parte[8].Inicie( 0,0,0, 0,0,0,       0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);
  parte[9].Inicie( 0,0,0, M_PI/2,0,0,  0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);
  parte[10].Inicie(0,0,0, 0,0,-M_PI/2, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);
  parte[11].Inicie(0,A,0, M_PI/2,0,0,  0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);
  parte[12].Inicie(0,A,0, 0,M_PI/2,0,  0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);
  parte[13].Inicie(A,0,0, 0,0,-M_PI/2, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);
  parte[14].Inicie(A,0,0, 0,0,0,       0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);
  parte[15].Inicie(0,0,A, M_PI/2,0,0,  0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);
  parte[16].Inicie(0,0,A, 0,0,-M_PI/2, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);
  parte[17].Inicie(A,A,0, 0,M_PI/2,0,  0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);
  parte[18].Inicie(0,A,A, M_PI/2,0,0,  0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);
  parte[19].Inicie(A,0,A, 0,0,-M_PI/2, 0,0,0,0,0,0,0,0,0, V,Vy,0,theta0,omega0,mc,Ic);

  //PLANOS
  parte[20].Inicie(0,0,0,0,0,0, -R,0,0,  -R,A,0,  -R,0,A,  V,Vy,0,theta0,omega0,mp,Ip);
  parte[21].Inicie(0,0,0,0,0,0, 0,0,-R,  0,A,-R,  A,0,-R,  V,Vy,0,theta0,omega0,mp,Ip);
  parte[22].Inicie(0,0,0,0,0,0, 0,-R,0,  A,-R,0,  0,-R,A,  V,Vy,0,theta0,omega0,mp,Ip);
  parte[23].Inicie(0,0,0,0,0,0, 0,0,A+R, A,0,A+R, 0,A,A+R, V,Vy,0,theta0,omega0,mp,Ip);
  parte[24].Inicie(0,0,0,0,0,0, A+R,0,0, A+R,0,A, A+R,A,0, V,Vy,0,theta0,omega0,mp,Ip);
  parte[25].Inicie(0,0,0,0,0,0, 0,A+R,0, A,A+R,0, 0,A+R,A, V,Vy,0,theta0,omega0,mp,Ip);

  parte[26].Inicie(0,0,0,0,0,0, -L/2,-L/2,-L/2, -L/2,L/2,-L/2, L/2,-L/2,-L/2, 0,0,0,0,0,mp_grande,Ip_grande);
  
int i;
  //Dibuja cubo y plano
  for (t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    
    if (tdibujo>tmax/Ndibujos){
      
      InicieCuadro();
 
      for(i=0;i<8;i++) parte[i].Dibujar_esfera();
      for(i=8;i<20;i++) parte[i].Dibujar_cilindro();
      for(i=20;i<N;i++) parte[i].Dibujar_plano();
      
      TermineCuadro();
      tdibujo=0;
      
    }
   

      //Muevase con Omelyan FR.
      for(i=0;i<N;i++){
	parte[i].Mueva_r(dt, Zeta);
      }

      Newton.CalculeTodasLasFuerzas(parte,t);
    
      for(i=0;i<N;i++){
        parte[i].Mueva_V(dt, (1-2*Lambda)/2);
      }
      for(i=0;i<N;i++){
        parte[i].Mueva_r(dt, Xi);
      }
   
      Newton.CalculeTodasLasFuerzas(parte,t);
    
      for(i=0;i<N;i++){
        parte[i].Mueva_V(dt, Lambda);
      }
      for(i=0;i<N;i++){
        parte[i].Mueva_r(dt, 1-2*(Xi+Zeta)); 
      }
      Newton.CalculeTodasLasFuerzas(parte,t);
    
      for(i=0;i<N;i++){
        parte[i].Mueva_V(dt, Lambda);
      }
      for(i=0;i<N;i++){
        parte[i].Mueva_r(dt, Xi);
      }
    
      Newton.CalculeTodasLasFuerzas(parte,t);
    
      for(i=0;i<N;i++){
        parte[i].Mueva_V(dt, (1-2*Lambda)/2);
      }
      for(i=0;i<N;i++){
        parte[i].Mueva_r(dt, Zeta);
      }

  }
  
  return 0;

}
