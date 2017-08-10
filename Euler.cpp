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


void UnPasoDeEuler(double &x, double &t, double dt){
  double dx=f(x, t)*dt;
  x+=dx;
}



int main(void){

  double x=1, t=0;
  double dt=0.01;

  for(t=0; t<10; t+=dt){
    cout<<t<<" "<<x<<" "<<fExacta(t)<<endl;
    UnPasoDeEuler(x, t, dt);
  }
    
  
  return 0;

}
