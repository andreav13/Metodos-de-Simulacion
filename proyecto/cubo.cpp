#include <iostream>
#include <fstream>
#include <cmath> 
#include "Vector.h"
#include "Random64.h"

using namespace std;

void Dibujar_esfera(double x, double y, double z, double R){
  cout<<", "<<x<<"+"<<R<<"*sin(u)*cos(v),"<<y<<"+"<<R<<"*sin(u)*sin(v),"<<z<<"+"<<R<<"*cos(u)";
}

void Dibujar_cilindro_rot_y_z_x(double x, double y, double z, double R, double phi, double alpha, double beta){

  cout<<", "<<x<<"+"<<R<<"*cos(v)*cos("<<phi<<")*cos("<<alpha<<")"<<"+t(u)*sin("<<phi<<")*cos("<<alpha<<")-"<<R<<"*sin(v)*sin("<<alpha<<"),"<<y<<"+"<<R<<"*cos(v)*cos("<<phi<<")*sin("<<alpha<<")*cos("<<beta<<")"<<"+t(u)*sin("<<phi<<")*sin("<<alpha<<")*cos("<<beta<<")+"<<R<<"*sin(v)*cos("<<alpha<<")*cos("<<beta<<")+"<<R<<"*cos(v)*sin("<<phi<<")*sin("<<beta<<")-cos("<<phi<<")*t(u)*sin("<<beta<<"),"<<z<<"+"<<R<<"*cos(v)*cos("<<phi<<")*sin("<<alpha<<")*sin("<<beta<<")+sin("<<phi<<")*t(u)*sin("<<alpha<<")*sin("<<beta<<")+"<<R<<"*sin(v)*cos("<<alpha<<")*sin("<<beta<<")-"<<R<<"*cos(v)*sin("<<phi<<")*cos("<<beta<<")+cos("<<phi<<")*t(u)*cos("<<beta<<")";
  
}

void Dibujar_plano(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3){
  cout<<", "<<a1<<"+m(u)*("<<b1<<"-"<<a1<<")+n(v)*("<<c1<<"-"<<a1<<"),"<<a2<<"+m(u)*("<<b2<<"-"<<a2<<")+n(v)*("<<c2<<"-"<<a2<<"),"<<a3<<"+m(u)*("<<b3<<"-"<<a3<<")+n(v)*("<<c3<<"-"<<a3<<")";
}

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output '1cubo.gif'"<<endl;
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
  double a=10, L=100., R=1, sqrt2=pow(2,0.5);

  InicieAnimacion();

  for (t=0;t<tmax;t+=dt){
      InicieCuadro();

      //DIBUJA CUBO
      
      Dibujar_esfera(0,  0,  0, R);
      Dibujar_esfera(a, 0, 0, R);
      Dibujar_esfera(0, a, 0, R);
      Dibujar_esfera(0, 0, a, R);
      Dibujar_esfera(0, a, a, R);
      Dibujar_esfera(a, 0, a, R);
      Dibujar_esfera(a, a, 0, R);
      Dibujar_esfera(a, a, a, R);

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
