#include <iostream>
#include <fstream>
#include <cmath> 
#include "Vector.h"
#include "Random64.h"

using namespace std;

void Dibujar_esfera(double x, double y, double z, double R){
  cout<<", "<<x<<"+"<<R<<"*sin(u)*cos(v),"<<y<<"+"<<R<<"*sin(u)*sin(v),"<<z<<"+"<<R<<"*cos(u)";
}

void Dibujar_cilindro_rot_y_z_x(double x, double y, double z, double R, double L, double phi, double alpha, double beta){

  string s;
  if(L==10) s="t";
  if(L==14) s="s";
  
  cout<<", "<<x<<"+"<<R<<"*cos(v)*cos("<<phi<<")*cos("<<alpha<<")"<<"+"<<s<<"(u)*sin("<<phi<<")*cos("<<alpha<<")-"<<R<<"*sin(v)*sin("<<alpha<<"),"<<y<<"+"<<R<<"*cos(v)*cos("<<phi<<")*sin("<<alpha<<")*cos("<<beta<<")"<<"+"<<s<<"(u)*sin("<<phi<<")*sin("<<alpha<<")*cos("<<beta<<")+"<<R<<"*sin(v)*cos("<<alpha<<")*cos("<<beta<<")+"<<R<<"*cos(v)*sin("<<phi<<")*sin("<<beta<<")-cos("<<phi<<")*"<<s<<"(u)*sin("<<beta<<"),"<<z<<"+"<<R<<"*cos(v)*cos("<<phi<<")*sin("<<alpha<<")*sin("<<beta<<")+sin("<<phi<<")*"<<s<<"(u)*sin("<<alpha<<")*sin("<<beta<<")+"<<R<<"*sin(v)*cos("<<alpha<<")*sin("<<beta<<")-"<<R<<"*cos(v)*sin("<<phi<<")*cos("<<beta<<")+cos("<<phi<<")*"<<s<<"(u)*cos("<<beta<<")";
  
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
}

void TermineCuadro(void){
  cout<<endl;
}

int main(void){

  double t, dt=1e-3, tmax=3*dt;
  double a=10, R=1, L1=10, L2=14;

  InicieAnimacion();

  for (t=0;t<tmax;t+=dt){
      InicieCuadro();

      Dibujar_esfera(0,  0,  0, R);
      Dibujar_esfera(a, 0, 0, R);
      Dibujar_esfera(0, a, 0, R);
      Dibujar_esfera(0, 0, a, R);

      Dibujar_cilindro_rot_y_z_x(0, 0, 0, R, L1, 0, 0, 0);
      Dibujar_cilindro_rot_y_z_x(0, 0, 0, R, L1, M_PI/2, 0, 0);
      Dibujar_cilindro_rot_y_z_x(0, 0, 0, R, L1, 0, 0, -M_PI/2);
      Dibujar_cilindro_rot_y_z_x(0, a, 0, R, L2, M_PI/2, -M_PI/4, 0);
      Dibujar_cilindro_rot_y_z_x(0, a, 0, R, L2, 0, -M_PI/4, M_PI/4);
      Dibujar_cilindro_rot_y_z_x(a, 0, 0, R, L2, -M_PI/4, 0, 0);
      
      TermineCuadro();
  }
  
  return 0;

}
