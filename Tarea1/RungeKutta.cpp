#include <iostream> //print
#include <fstream> //save to file
#include <cmath> //math functions

using namespace std;


double fExacta(double t){
  return exp(t);
}


double f(double x, double t){
  return x;
}


void UnPasoDeRungeKutta4(double &x, double &t, double dt){
  double dx1, dx2, dx3, dx4;
  dx1=dt*f(x, t);
  dx2=dt*f(x+dx1/2., t+dt/2.);
  dx3=dt*f(x+dx2/2., t+dt/2.);
  dx4=dt*f(x+dx3, t+dt);
  
  x += dx1/6. + dx2/3. + dx3/3. + dx4/6.;
}



int main(void){

  double x=1, t=0;
  double dt=0.1;

  for(t=0; t<10; t+=dt){
    cout<<t<<" "<<x<<" "<<fExacta(t)<<endl;
    UnPasoDeRungeKutta4(x, t, dt);
  }
    
  
  return 0;

}
