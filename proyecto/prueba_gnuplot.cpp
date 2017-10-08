#include <iostream>
#include <fstream> //save to file
#include <cmath> //math functions
#include "Vector.h"
#include "Random64.h"

using namespace std;

void Dibujar_circulo(double x, double y, double z, double R){
  //cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
  //  <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";

  cout<<", "<<x<<"+"<<R<<"*sin(u)*cos(v),"<<y<<"+"<<R<<"*sin(u)*sin(v),"<<z<<"+"<<R<<"*cos(u)";
}

void Dibujar_cilindro(double x, double y, double z, double R){

  cout<<", "<<x<<"+"<<R<<"*cos(v),"<<y<<"+"<<R<<"*sin(v),"<<z<<"+"<<"t(u)";
}

void Dibujar_cilindro_rot_z(double x, double y, double z, double R, double phi){

  cout<<", "<<x<<"+"<<R<<"*cos(v)*cos("<<phi<<")"<<"-"<<R<<"*sin(v)*sin("<<phi<<"),"<<y<<"+"<<R<<"*sin(v)*cos("<<phi<<")"<<"+"<<R<<"*cos(v)*sin("<<phi<<"),"<<z<<"+"<<"t(u)";
}

void Dibujar_cilindro_rot_y(double x, double y, double z, double R, double phi){

  cout<<", "<<x<<"+"<<R<<"*cos(v)*cos("<<phi<<")"<<"+"<<"t(u)*sin("<<phi<<"),"<<y<<"+"<<R<<"*sin(v),"<<z<<"-"<<R<<"*cos(v)*sin("<<phi<<")"<<"+"<<"t(u)*cos("<<phi<<")";
}

void Dibujar_cilindro_rot_x(double x, double y, double z, double R, double phi){

  cout<<", "<<x<<"+"<<R<<"*cos(v),"<<y<<"+"<<R<<"*sin(v)*cos("<<phi<<")"<<"-"<<"t(u)*sin("<<phi<<"),"<<z<<"+"<<R<<"*sin(v)*sin("<<phi<<")"<<"+"<<"t(u)*cos("<<phi<<")";
}

void Dibujar_cilindro_rot_y_z(double x, double y, double z, double R, double phi, double alpha){

  cout<<", "<<x<<"+"<<R<<"*cos(v)*cos("<<phi<<")*cos("<<alpha<<")"<<"+"<<"s(u)*sin("<<phi<<")*cos("<<alpha<<")-"<<R<<"*sin(v)*sin("<<alpha<<"),"<<y<<"+"<<R<<"*cos(v)*cos("<<phi<<")*sin("<<alpha<<")"<<"+"<<"s(u)*sin("<<phi<<")*sin("<<alpha<<")+"<<R<<"*sin(v)*cos("<<alpha<<"),"<<z<<"-"<<R<<"*cos(v)*sin("<<phi<<")"<<"+"<<"s(u)*cos("<<phi<<")";
}

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'prueba1.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-2:12]"<<endl;
  cout<<"set yrange [-2:12]"<<endl;
  cout<<"set zrange [-2:12]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set view equal xyz"<<endl;
  cout<<"set hidden3d front"<<endl;
  cout<<"set pm3d depthorder interpolate 0,0"<<endl;
  cout<<"set palette rgb 3,3,3"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set urange [0:pi]"<<endl;
  cout<<"set vrange [0:2*pi]"<<endl;
  cout<<"t(u)=10*u/pi"<<endl;
  cout<<"s(u)=14.14213562*u/pi"<<endl;
  cout<<"set isosamples 20"<<endl;
}

void InicieCuadro(void){
  cout<<"splot 0,0,0 ";
  /*cout<<" , "<<100.1/7<<"*t,0";
  cout<<" , "<<100.1/7<<"*t,100";
  cout<<" , 0,"<<100.1/7<<"*t";
  cout<<" , 100,"<<100.1/7<<"*t";*/
}

void TermineCuadro(void){
  cout<<endl;
}

int main(void){

  double t, dt=1e-3, tmax=3*dt;
  double a=10, R=1, L=sqrt(2*20*20);

  InicieAnimacion();

  for (t=0;t<tmax;t+=dt){
      InicieCuadro();
      
      /*Dibujar_circulo(a,  a,  a, R);
      Dibujar_circulo(a, -a, -a, R);
      Dibujar_circulo(-a, a, -a, R);
      Dibujar_circulo(-a, -a, a, R);*/

      Dibujar_circulo(0,  0,  0, R);
      Dibujar_circulo(a, 0, 0, R);
      Dibujar_circulo(0, a, 0, R);
      Dibujar_circulo(0, 0, a, R);

      //Dibujar_cilindro(0, 0, 0, R);
      Dibujar_cilindro_rot_z(0, 0, 0, R, M_PI/3);
      Dibujar_cilindro_rot_y(0, 0, 0, R, M_PI/2);
      Dibujar_cilindro_rot_x(0, 0, 0, R, -M_PI/2);
      Dibujar_cilindro_rot_y_z(0, a, 0, R, M_PI/2, -M_PI/4);
      
      TermineCuadro();
  }
  
  return 0;

}
