#include <iostream>
#include <fstream> 
#include <cmath> 
#include "Random64.h"

using namespace std;

const double Deltat=0.01, L=100, T=300, MASA= 22.8916, e=1, D=1.132, kb=0.826;//unidades: picosegundos, Amstromg, Kelvin, u.m.a.
const double Gamma=kb*T/(MASA*D), sigma=sqrt(2*D*Deltat), dtU2mGamma=Deltat/(2*MASA*Gamma);

const int N=4000;

class Cuerpo;

//Clase Cuerpo
class Cuerpo{
private:
  double x,xold,Fx,Fxold,m,R;
  
public:
  void Inicie(double x0, double m0, double R0);
  void CalculeFuerza(double E);
  void Dibujese(void);
  void Muevase(Crandom & ran64);
  double Getx(void){return x;};
  double Getj(void){return (x-xold)/Deltat*e;};
};


//Funciones de la clase cuerpo
void Cuerpo::Inicie(double x0, double m0, double R0){
  x=xold=x0; m=m0; R=R0; Fx=0;
}


void Cuerpo::CalculeFuerza(double E){
  Fxold=Fx; Fx=e*E;
}


void Cuerpo::Muevase(Crandom & ran64){
  xold=x;
  x+=dtU2mGamma*(3*Fx-Fxold) + ran64.gauss(0,sigma); 
}

void Cuerpo::Dibujese(void){
  cout<<", "<<x<<"+"<<R<<"*cos(t),"<<0<<"+"<<R<<"*sin(t)";
}


//Funciones globales
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'DinBrownianaMuchos.gif'"<<endl;
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

double GetSigma2(Cuerpo * Na){
  double xprom, sigma2;
  int i;
  
  //calculo xprom
  for(xprom=0,i=0;i<N;i++){
    xprom+=Na[i].Getx();
  }  xprom/=N;
  
  //calculo sigma
  for(sigma2=0,i=0;i<N;i++){
    sigma2+=pow(Na[i].Getx()-xprom,2);
  }sigma2/=(N-1);
  
  return sigma2;

}

double Getjprom(Cuerpo * Na){
  double jprom; int i;
  for(jprom=0,i=0;i<N;i++){
    jprom+=Na[i].Getj();
  }jprom/=N;
  return jprom;
}


//Main
int main(void){

  double t,  tmax=40;
  Crandom ran64(1);
  double E=200, sigmax0=5;
  Cuerpo Na[N];   //Sodio
  int i;

  
  //InicieAnimacion();
  
  //(x0, m0, R0)
  for(i=0;i<N;i++){Na[i].Inicie(ran64.gauss(L/2, sigmax0), MASA, 4);}
  for(i=0;i<N;i++){Na[i].CalculeFuerza(E);}

  for(t=0;t<tmax;t+=Deltat){
    cout<<t<<" "<<Getjprom(Na)<<endl;
      
    //InicieCuadro();
    //Na[0].Dibujese();
    //TermineCuadro();
     
    for(i=0;i<N;i++){Na[i].CalculeFuerza(E);}
    for(i=0;i<N;i++){Na[i].Muevase(ran64);}
  }

    //IMPRIMA
    //for(i=0;i<N;i++){cout<<Na[i].Getx()<<endl;}
  
  return 0;

}
