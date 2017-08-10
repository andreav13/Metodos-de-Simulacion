#include <iostream>
#include <fstream> //save to file
#include <cmath> //math functions

using namespace std;

const double GM=1;
class Cuerpo;

//Clase Cuerpo
class Cuerpo{
private:
  double x,y,Vx,Vy,Fx,Fy,m,R;
  
public:
  void Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0);
  void CalculeFuerza(void);
  void Muevase(double dt);
  double Getx(void){return x;};
  double Gety(void){return y;};
};


//Funciones de la clase cuerpo
void Cuerpo::Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0){
  x=x0; y=y0; Vx=Vx0; Vy=Vy0; m=m0; R=R0;
}


void Cuerpo::CalculeFuerza(void){
  double r2=x*x+y*y;
  double aux=GM*m*pow(r2,-1.5);
  Fx=-aux*x; Fy=-aux*y;
}


void Cuerpo::Muevase(double dt){
  x+=Vx*dt; y+=Vy*dt;
  Vx+=Fx*dt/m; Vy+=Fy*dt/m;
}


//Main
int main(void){

  double t, dt=0.01;
  Cuerpo Planeta;

  //(x0, y0, Vx0, Vy0, m0, R0)
  Planeta.Inicie(0, 0, 12, 16, 0.457, 0.15);
  
  for(t=0;t<3.5;t+=dt){
    cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
    Planeta.CalculeFuerza();
    Planeta.Muevase(dt);
  }
  
  return 0;

}
