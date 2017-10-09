#include <iostream>
#include <fstream>
#include <cmath> 
#include "Vector.h"
#include "Random64.h"

using namespace std;

const double A=10., L=100., R=1.;
const double K=1e4, Gamma=50, Kcundall=10, MU=0.4;
const double g=9.8;

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;
  
class Cuerpo;
class Colisionador;

//Clase Cuerpo
class Cuerpo{
private:
  vector3D r,rot,a,b,c,V,F,omega,tau;
  double m,R,theta,I;
  
public:
  void Inicie(double x0, double y0, double z0, double phi0, double alpha0, double beta0, double a10, double a20, double a30, double b10, double b20, double b30,  double c10, double c20, double c30, double Vx0, double Vy0, double Vz0,double theta0, double omega0, double m0, double R0);
  void BorreFuerzayTorque(void);
  void AgregueFuerza(vector3D F0);
  void AgregueTorque(vector3D tau0);
  void Mueva_r(double dt, double Constante);
  void Mueva_V(double dt, double Constante);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  double Getz(void){return r.z();};
  void Dibujar_esfera(void);
  void Dibujar_cilindro(void);
  void Dibujar_plano(void);
  friend class Colisionador;
};


//Clase Colisionador
class Colisionador{
private:
  vector3D ele[N+4][N+4];
  bool EstoyEnColision[N+4][N+4];
public:
  void Inicie(void);
  void CalculeTodasLasFuerzas(Cuerpo* Grano, double dt);
  void CalculeLaFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2, vector3D & ele, bool & EstoyEnColision, double dt);
};
  
//Funciones de la clase cuerpo
void Cuerpo::Inicie(double x0, double y0, double z0, double phi0, double alpha0, double beta0, double a10, double a20, double a30, double b10, double b20, double b30,  double c10, double c20, double c30, double Vx0, double Vy0, double Vz0,double theta0, double omega0, double m0, double R0){
  r.cargue(x0,y0,z0);
  rot.cargue(phi0,alpha0,beta0);
  a.cargue(a10,a20,a30);
  b.cargue(b10,b20,b3);
  c.cargue(c10,c20,c30);
  V.cargue(Vx0,Vy0,Vz0);
  theta=theta0;
  omega.cargue(0,0,omega0);
  m=m0; R=R0; I=m*a*a/6.;
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


void Cuerpo::Dibujar_esfera(void){
  cout<<", "<<r.x()<<"+"<<R<<"*sin(u)*cos(v),"<<r.y()<<"+"<<R<<"*sin(u)*sin(v),"<<r.z()<<"+"<<R<<"*cos(u)";
}

void Cuerpo::Dibujar_cilindro(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(v)*cos("<<rot.x()<<")*cos("<<rot.y()<<")"<<"+t(u)*sin("<<rot.x()<<")*cos("<<rot.y()<<")-"<<R<<"*sin(v)*sin("<<rot.y()<<"),"<<y<<"+"<<R<<"*cos(v)*cos("<<rot.x()<<")*sin("<<rot.y()<<")*cos("<<rot.z()<<")"<<"+t(u)*sin("<<rot.x()<<")*sin("<<rot.y()<<")*cos("<<rot.z()<<")+"<<R<<"*sin(v)*cos("<<rot.y()<<")*cos("<<rot.z()<<")+"<<R<<"*cos(v)*sin("<<rot.x()<<")*sin("<<rot.z()<<")-cos("<<rot.x()<<")*t(u)*sin("<<rot.z()<<"),"<<z<<"+"<<R<<"*cos(v)*cos("<<rot.x()<<")*sin("<<rot.y()<<")*sin("<<rot.z()<<")+sin("<<rot.x()<<")*t(u)*sin("<<rot.y()<<")*sin("<<rot.z()<<")+"<<R<<"*sin(v)*cos("<<rot.y()<<")*sin("<<rot.z()<<")-"<<R<<"*cos(v)*sin("<<rot.x()<<")*cos("<<rot.z()<<")+cos("<<rot.x()<<")*t(u)*cos("<<rot.z()<<")";
}

void Cuerpo::Dibujar_plano(void){
  cout<<", "<<a.x()<<"+m(u)*("<<b.x()<<"-"<<a.x()<<")+n(v)*("<<c.x()<<"-"<<a.x()<<"),"<<a.y()<<"+m(u)*("<<b.y()<<"-"<<a.y()<<")+n(v)*("<<c.y()<<"-"<<a.y()<<"),"<<a.z()<<"+m(u)*("<<b.z()<<"-"<<a.z()<<")+n(v)*("<<c.z()<<"-"<<a.z()<<")";
}


//Funciones de la clase colisionador

void Colisionador::Inicie(void){
  int i,j;
  for(i=0;i<N;i++){
    for(j=i+1;j<N+4;j++){
      ele[i][j].cargue(0,0,0);
      EstoyEnColision[i][j]=false;
    }
  }
}

void Colisionador::CalculeTodasLasFuerzas(Cuerpo* Grano, double dt){
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
      CalculeLaFuerzaEntre(Grano[i], Grano[j], ele[i][j], EstoyEnColision[i][j], dt);
    }
  }
}


void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2, vector3D & ele, bool & EstoyEnColision, double dt){
  vector3D F2,Fn,Ft,runitario,n,Vc,Vcn,Vct,t,r21 = Grano2.r-Grano1.r;
  double s,m1,m2,R1,R2,m12,componenteVcn,normaVct,componenteFn,normaFt,Ftmax;
  double ERFF=1e-8,d21=norma(r21);
  s=Grano1.R+Grano2.R-d21;

  if(s>0){    //SI SE CHOCAN

    //Geometria y dinamica del contacto
    m1=Grano1.m; m2=Grano2.m; m12=(m1*m2)/(m1+m2);
    R1=Grano1.R; R2=Grano2.R;
    n=r21/d21;

    //Calcular velocidad de contacto y el vector tangente
    Vc=(Grano2.V-Grano1.V)-(Grano2.omega^n)*R2-(Grano1.omega^n)*R1;  // ^ = prod cruz
    //velocidad del punto de contacto
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

    //FUERZAS TANGENCIALES
    //Fuerza estática
    ele+=(Vct*dt);
    Ft=ele*(-Kcundall);
    //Fuerza cinética
    Ftmax=MU*componenteFn; normaFt=norma(Ft);
    if(normaFt>Ftmax) Ft=ele*(-Ftmax/norma(ele));
    
    F2=Fn+Ft;
    Grano1.AgregueFuerza(F2*(-1)); Grano2.AgregueFuerza(F2);
    Grano1.AgregueTorque((n*(-R2))^Ft); Grano2.AgregueTorque((n*R1)^(Ft*(-1)));

    EstoyEnColision=true;
  }
  
  else if(EstoyEnColision==true){
    ele.cargue(0,0,0); EstoyEnColision=false;
  }
  
}



void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output '2cubo.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-60:60]"<<endl;
  cout<<"set yrange [-60:60]"<<endl;
  cout<<"set zrange [-60:12]"<<endl;
  cout<<"set term gif size 900,700"<<endl;
  cout<<"set size 1.1,1.1"<<endl;
  cout<<"set size ratio 1"<<endl;
  cout<<"set view equal xyz"<<endl;
  cout<<"set hidden3d front"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set urange [0:pi]"<<endl;
  cout<<"set vrange [0:2*pi]"<<endl;
  cout<<"t(u)=10*u/pi"<<endl;
  cout<<"m(u)=u/pi"<<endl;
  cout<<"n(v)=v/(2*pi)"<<endl;
  cout<<"set isosamples 8,8"<<endl;
}

void InicieCuadro(void){
  cout<<"splot 0,0,0 ";
}

void TermineCuadro(void){
  cout<<endl;
}

int main(void){

  double t, dt=1e-3, tmax=3*dt;
  Cuerpo parte[26];
  Colisionador Choque;
  Crandom ran64(1);

  double m0=1, V=10, theta0=M_PI, omega0=10;
  
  InicieAnimacion();

  //x0, y0, z0,   phi0,   alpha0,   beta0,   a10,   a20,   a30,   b10,   b20,   b30,    c10,   c20,   c30,   Vx0,   Vy0,   Vz0,  theta0,   omega0,   m0,   R0

  //ESFERAS
  parte[0].Inicie(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,V,0,0,theta0,omega0,m0,R);
  parte[1].Inicie(A,0,0,0,0,0,0,0,0,0,0,0,0,0,0,V,0,0,theta0,omega0,m0,R);
  parte[2].Inicie(0,A,0,0,0,0,0,0,0,0,0,0,0,0,0,V,0,0,theta0,omega0,m0,R);
  parte[3].Inicie(0,0,A,0,0,0,0,0,0,0,0,0,0,0,0,V,0,0,theta0,omega0,m0,R);
  parte[4].Inicie(0,A,A,0,0,0,0,0,0,0,0,0,0,0,0,V,0,0,theta0,omega0,m0,R);
  parte[5].Inicie(A,0,A,0,0,0,0,0,0,0,0,0,0,0,0,V,0,0,theta0,omega0,m0,R);
  parte[6].Inicie(A,A,0,0,0,0,0,0,0,0,0,0,0,0,0,V,0,0,theta0,omega0,m0,R);
  parte[7].Inicie(A,A,A,0,0,0,0,0,0,0,0,0,0,0,0,V,0,0,theta0,omega0,m0,R);

  //CILINDROS
  parte[8].Inicie();
  
  
  for (t=0;t<tmax;t+=dt){
      InicieCuadro();
    



      Dibujar_cilindro_rot_y_z_x(0, 0, 0, R, 0, 0, 0);
      Dibujar_cilindro_rot_y_z_x(0, 0, 0, R, M_PI/2, 0, 0);
      Dibujar_cilindro_rot_y_z_x(0, 0, 0, R, 0, 0, -M_PI/2);
      Dibujar_cilindro_rot_y_z_x(0, a, 0, R, M_PI/2, 0, 0);
      Dibujar_cilindro_rot_y_z_x(0, a, 0, R, 0, M_PI/2, 0);
      Dibujar_cilindro_rot_y_z_x(a, 0, 0, R, 0, 0, -M_PI/2);
      
      Dibujar_cilindro_rot_y_z_x(a, 0, 0, R, 0, 0, 0);
      Dibujar_cilindro_rot_y_z_x(0, 0, a, R, M_PI/2, 0, 0);
      Dibujar_cilindro_rot_y_z_x(0, 0, a, R, 0, 0, -M_PI/2);
      Dibujar_cilindro_rot_y_z_x(a, a, 0, R, 0, M_PI/2, 0);
      Dibujar_cilindro_rot_y_z_x(0, a, a, R, M_PI/2, 0, 0);
      Dibujar_cilindro_rot_y_z_x(a, 0, a, R, 0, 0, -M_PI/2);

      Dibujar_plano(-R,0,0, -R,a,0, -R,0,a);
      Dibujar_plano(0,0,-R, 0,a,-R, a,0,-R);
      Dibujar_plano(0,-R,0, a,-R,0, 0,-R,a);
      Dibujar_plano(0,0,a+R, a,0,a+R, 0,a,a+R);
      Dibujar_plano(a+R,0,0, a+R,0,a, a+R,a,0);
      Dibujar_plano(0,a+R,0, a,a+R,0, 0,a+R,a);

      //DIBUJA PLANO

      Dibujar_plano(-L/2,-L/2,-L/2, -L/2,L/2,-L/2, L/2,-L/2,-L/2);
      
      TermineCuadro();
  }
  
  return 0;

}
