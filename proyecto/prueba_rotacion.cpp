#include <iostream>
#include <fstream>
#include <cmath> 
#include "Vector.h"
#include "Random64.h"

using namespace std;

const double A=10., L=100., R=1.;
const double K=1e4, Gamma=50, Kcundall=10;
const double g=9.8;
const int N=27;

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;
  
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
  void Dibujar_esfera(void);
  void Dibujar_cilindro(void);
  void Dibujar_plano(void);
  void Rotar_esfera(double theta, double phi, double h, double k, double l);
  void Rotar_cilindro1(double theta, double phi, double h, double k, double l);
  void Rotar_cilindro2(double theta, double phi, double h, double k, double l);
  void Rotar_plano(double theta, double phi, double h, double k, double l);
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
  void CalculeLaFuerzaEntrePlanos(Cuerpo * Grano, vector3D & ele, bool & EstoyEnColision, double dt);
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
  omega+=tau*Constante*dt/I;
}


void Cuerpo::Dibujar_esfera(void){
  cout<<", "<<r.x()<<"+"<<R<<"*sin(u)*cos(v),"<<r.y()<<"+"<<R<<"*sin(u)*sin(v),"<<r.z()<<"+"<<R<<"*cos(u)";
}

void Cuerpo::Dibujar_cilindro(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(v)*cos("<<rot.x()<<")*cos("<<rot.y()<<")"<<"+t(u)*sin("<<rot.x()<<")*cos("<<rot.y()<<")-"<<R<<"*sin(v)*sin("<<rot.y()<<"),"<<r.y()<<"+"<<R<<"*cos(v)*cos("<<rot.x()<<")*sin("<<rot.y()<<")*cos("<<rot.z()<<")"<<"+t(u)*sin("<<rot.x()<<")*sin("<<rot.y()<<")*cos("<<rot.z()<<")+"<<R<<"*sin(v)*cos("<<rot.y()<<")*cos("<<rot.z()<<")+"<<R<<"*cos(v)*sin("<<rot.x()<<")*sin("<<rot.z()<<")-cos("<<rot.x()<<")*t(u)*sin("<<rot.z()<<"),"<<r.z()<<"+"<<R<<"*cos(v)*cos("<<rot.x()<<")*sin("<<rot.y()<<")*sin("<<rot.z()<<")+sin("<<rot.x()<<")*t(u)*sin("<<rot.y()<<")*sin("<<rot.z()<<")+"<<R<<"*sin(v)*cos("<<rot.y()<<")*sin("<<rot.z()<<")-"<<R<<"*cos(v)*sin("<<rot.x()<<")*cos("<<rot.z()<<")+cos("<<rot.x()<<")*t(u)*cos("<<rot.z()<<")";
}

void Cuerpo::Dibujar_plano(void){
  cout<<", "<<a.x()<<"+m(u)*("<<b.x()<<"-"<<a.x()<<")+n(v)*("<<c.x()<<"-"<<a.x()<<"),"<<a.y()<<"+m(u)*("<<b.y()<<"-"<<a.y()<<")+n(v)*("<<c.y()<<"-"<<a.y()<<"),"<<a.z()<<"+m(u)*("<<b.z()<<"-"<<a.z()<<")+n(v)*("<<c.z()<<"-"<<a.z()<<")";
}


void Cuerpo::Rotar_esfera(double theta, double phi, double h, double k, double l){
  double r1,r2,r3;

  r1=(r.x()-h)*cos(theta)-(r.y()-k)*sin(theta)+h;
  r2=(r.x()-h)*sin(theta)+(r.y()-k)*cos(theta)+k;
  r3=r.z();
  r.cargue(r1,r2,r3);

  r1=r.x();
  r2=(r.y()-k)*cos(phi)-(r.z()-l)*sin(phi)+k;
  r3=(r.y()-k)*sin(phi)+(r.z()-l)*cos(phi)+l;
  r.cargue(r1,r2,r3);

}

void Cuerpo::Rotar_cilindro1(double theta, double phi, double h, double k, double l){
  double r1,r2,r3;

  r1=(r.x()-h)*cos(theta)-(r.y()-k)*sin(theta)+h;
  r2=(r.x()-h)*sin(theta)+(r.y()-k)*cos(theta)+k;
  r3=r.z();
  r.cargue(r1,r2,r3);

  /*r1=(r.x()-h)*cos(phi)+(r.z()-l)*sin(phi)+h;
  r2=r.y();
  r3=-(r.x()-h)*sin(phi)+(r.z()-l)*cos(phi)+l;
  r.cargue(r1,r2,r3);*/

  r1=r.x();
  r2=(r.y()-k)*cos(phi)-(r.z()-l)*sin(phi)+k;
  r3=(r.y()-k)*sin(phi)+(r.z()-l)*cos(phi)+l;
  r.cargue(r1,r2,r3);

  rot.cargue(rot.x(), rot.y()+theta, rot.z()+phi);

}

void Cuerpo::Rotar_cilindro2(double theta, double phi, double h, double k, double l){
  double r1,r2,r3;

  r1=(r.x()-h)*cos(theta)-(r.y()-k)*sin(theta)+h;
  r2=(r.x()-h)*sin(theta)+(r.y()-k)*cos(theta)+k;
  r3=r.z();
  r.cargue(r1,r2,r3);

  /*r1=(r.x()-h)*cos(phi)+(r.z()-l)*sin(phi)+h;
  r2=r.y();
  r3=-(r.x()-h)*sin(phi)+(r.z()-l)*cos(phi)+l;
  r.cargue(r1,r2,r3);*/

  r1=r.x();
  r2=(r.y()-k)*cos(phi)-(r.z()-l)*sin(phi)+k;
  r3=(r.y()-k)*sin(phi)+(r.z()-l)*cos(phi)+l;
  r.cargue(r1,r2,r3);
  
  
  rot.cargue(rot.x()-theta, rot.y(), rot.z()+phi);
}

void Cuerpo::Rotar_plano(double theta, double phi, double h, double k, double l){
  double a1,a2,a3,b1,b2,b3,c1,c2,c3;
  
  a1=(a.x()-h)*cos(theta)-(a.y()-k)*sin(theta)+h;
  a2=(a.x()-h)*sin(theta)+(a.y()-k)*cos(theta)+k;
  a3=a.z();
  a.cargue(a1,a2,a3);

  b1=(b.x()-h)*cos(theta)-(b.y()-k)*sin(theta)+h;
  b2=(b.x()-h)*sin(theta)+(b.y()-k)*cos(theta)+k;
  b3=b.z();
  b.cargue(b1,b2,b3);

  c1=(c.x()-h)*cos(theta)-(c.y()-k)*sin(theta)+h;
  c2=(c.x()-h)*sin(theta)+(c.y()-k)*cos(theta)+k;
  c3=c.z();
  c.cargue(c1,c2,c3);


  a1=a.x();
  a2=(a.y()-k)*cos(phi)-(a.z()-l)*sin(phi)+k;
  a3=(a.y()-k)*sin(phi)+(a.z()-l)*cos(phi)+l;
  a.cargue(a1,a2,a3);

  b1=b.x();
  b2=(b.y()-k)*cos(phi)-(b.z()-l)*sin(phi)+k;
  b3=(b.y()-k)*sin(phi)+(b.z()-l)*cos(phi)+l;
  b.cargue(b1,b2,b3);

  c1=c.x();
  c2=(c.y()-k)*cos(phi)-(c.z()-l)*sin(phi)+k;
  c3=(c.y()-k)*sin(phi)+(c.z()-l)*cos(phi)+l;
  c.cargue(c1,c2,c3);

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
  int i;
  for(i=0;i<N;i++){
    Grano[i].BorreFuerzayTorque();
  }
  
  //Agregue la fuerza de la gravedad
  vector3D g_vector; g_vector.cargue(0,0,g);
  for(i=0;i<N-1;i++){
    Grano[i].AgregueFuerza(-(Grano[i].m)*g_vector);
  }
  
  //Calcular todas las fuerzas entre parejas de partes
  /*for(i=0;i<N;i++){
      CalculeLaFuerzaEntre(Grano[i], Grano[26], ele[i][26], EstoyEnColision[i][26], dt);
    }*/
  CalculeLaFuerzaEntrePlanos(Grano, ele[21][26], EstoyEnColision[21][26], dt);
}

/*
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
  
}*/

void Colisionador::CalculeLaFuerzaEntrePlanos(Cuerpo * Grano, vector3D & ele, bool & EstoyEnColision, double dt){
  vector3D F2,Fn,runitario,n,Vc,Vcn;
  double m1,m2,m12,componenteVcn,componenteFn;
  double ERFF=1e-8, s=-Grano[21].Getaz()-L/2;
  
  if(Grano[21].Getaz()<-L/2 and Grano[21].Getbz()<-L/2 and Grano[21].Getcz()<-L/2){    //SI SE CHOCAN
    //Geometria y dinamica del contacto
    n.cargue(0,0,-1);

    //Calcular velocidad de contacto 
    Vc=Grano[21].V; 
    componenteVcn=Vc*n;Vcn=n*componenteVcn;
    
    //FUERZAS NORMALES
    //Fuerza de Hertz
    componenteFn=K*pow(s,1.5);
    //Disipación Plástica
    //componenteFn-=m12*sqrt(s)*Gamma*componenteVcn;
    if(componenteFn<0) componenteFn=0;
    Fn=n*componenteFn;
    
    F2=Fn; int i;
    for(i=0;i<N-1;i++){
    Grano[i].AgregueFuerza(F2*(-1));
    }

    EstoyEnColision=true;
  }
  
}



void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'rot_cubo.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-60:60]"<<endl;
  cout<<"set yrange [-60:60]"<<endl;
  cout<<"set zrange [-60:12]"<<endl;
  //cout<<"set term gif size 900,700"<<endl;
  cout<<"set size 1.1,1.1"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set view equal xyz"<<endl;
  cout<<"set hidden3d front"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set urange [0:pi]"<<endl;
  cout<<"set vrange [0:2*pi]"<<endl;
  cout<<"t(u)=10*u/pi"<<endl;
  cout<<"m(u)=u/pi"<<endl;
  cout<<"n(v)=v/(2*pi)"<<endl;
  cout<<"set isosamples 8"<<endl;
}

void InicieCuadro(void){
  cout<<"splot 0,0,0";
}

void TermineCuadro(void){
  cout<<endl;
}

int main(void){

  double t, dt=1e-2, tmax=12;
  int Ndibujos;
  double tdibujo;
  Cuerpo parte[N];
  Colisionador Newton;
  Crandom ran64(1);

  double me=1,mc=1,mp=1,mp_grande=1000;
  double Ie=2/5.*me*R*R, Ic=1/2.*mc*R*R, Ip=1/6.*mp*A*A, Ip_grande=1/6.*mp_grande*L*L;
  double V=0, Vy=0, theta0=M_PI, omega0=10;

  Ndibujos=300;
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

  /*parte[20].Inicie(0,0,0,0,0,0, -0,0,0,  -0,A,0,  -0,0,A,  V,Vy,0,theta0,omega0,mp,Ip);
  parte[21].Inicie(0,0,0,0,0,0, 0,0,-0,  0,A,-0,  A,0,-0,  V,Vy,0,theta0,omega0,mp,Ip);
  parte[22].Inicie(0,0,0,0,0,0, 0,-0,0,  A,-0,0,  0,-0,A,  V,Vy,0,theta0,omega0,mp,Ip);
  parte[23].Inicie(0,0,0,0,0,0, 0,0,A+0, A,0,A+0, 0,A,A+0, V,Vy,0,theta0,omega0,mp,Ip);
  parte[24].Inicie(0,0,0,0,0,0, A+0,0,0, A+0,0,A, A+0,A,0, V,Vy,0,theta0,omega0,mp,Ip);
  parte[25].Inicie(0,0,0,0,0,0, 0,A+0,0, A,A+0,0, 0,A+0,A, V,Vy,0,theta0,omega0,mp,Ip);*/

  parte[26].Inicie(0,0,0,0,0,0, -L/2,-L/2,-L/2, -L/2,L/2,-L/2, L/2,-L/2,-L/2, 0,0,0,0,0,mp_grande,Ip_grande);


  int i;

  //Dibuja cubo y plano
  for (t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    if (tdibujo>tmax/Ndibujos){
      
      InicieCuadro();
 
      for(i=0;i<8;i++)  parte[i].Dibujar_esfera();
      for(i=8;i<20;i++) parte[i].Dibujar_cilindro();
      for(i=20;i<N;i++) parte[i].Dibujar_plano();

      TermineCuadro();
      tdibujo=0;
       
      }

    //Rota cubo
    Cuerpo A=parte[0], F=parte[7];
    double h=(A.Getx()+F.Getx())/2.,     k=(A.Gety()+F.Gety())/2.,     l=(A.Getz()+F.Getz())/2.;

    for (i=0;i<8;i++){parte[i].Rotar_esfera(M_PI/20,0,h,k,l);}
    for (i=8;i<20;i++)if(i%3!=1){parte[i].Rotar_cilindro1(M_PI/20,0,h,k,l);}
    for (i=8;i<20;i++)if(i%3==1){parte[i].Rotar_cilindro2(M_PI/20,0,h,k,l);}
    for (i=20;i<26;i++){parte[i].Rotar_plano(M_PI/20,0,h,k,l);}

    
    //Muevase con Omelyan FR.
    for(i=0;i<N;i++){
      parte[i].Mueva_r(dt, Zeta);
    }
    
    Newton.CalculeTodasLasFuerzas(parte, dt);
    
    for(i=0;i<N;i++){
      parte[i].Mueva_V(dt, (1-2*Lambda)/2);
    }
    for(i=0;i<N;i++){
      parte[i].Mueva_r(dt, Xi);
    }
    
    Newton.CalculeTodasLasFuerzas(parte, dt);
    
    for(i=0;i<N;i++){
      parte[i].Mueva_V(dt, Lambda);
    }
    for(i=0;i<N;i++){
      parte[i].Mueva_r(dt, 1-2*(Xi+Zeta)); 
    }
    Newton.CalculeTodasLasFuerzas(parte, dt);
    
    for(i=0;i<N;i++){
      parte[i].Mueva_V(dt, Lambda);
    }
    for(i=0;i<N;i++){
      parte[i].Mueva_r(dt, Xi);
    }
    
    Newton.CalculeTodasLasFuerzas(parte, dt);
    
    for(i=0;i<N;i++){
      parte[i].Mueva_V(dt, (1-2*Lambda)/2);
    }
    for(i=0;i<N;i++){
      parte[i].Mueva_r(dt, Zeta);
    }
    
  }
  
  return 0;
  
}
