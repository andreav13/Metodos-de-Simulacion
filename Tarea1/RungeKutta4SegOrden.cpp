#include <iostream> //print
#include <fstream> //save to file
#include <cmath> //math functions

using namespace std;

const double Omega0=1;

double fExacta(double t){
  return sin(Omega0*t); 
}


double f1(double x1, double x2, double t){
  return x2;
}


double f2(double x1, double x2, double t){
  return -Omega0*Omega0*x1;
}


void UnPasoDeRungeKutta4_2(double &x1, double &x2, double &t, double dt){
  double dx11, dx21, dx31, dx41;
  double dx12, dx22, dx32, dx42;
  
  dx11=dt*f1(x1, x2, t);
  dx12=dt*f2(x1, x2, t);
  
  dx21=dt*f1(x1+dx11/2., x2+dx12/2., t+dt/2.);
  dx22=dt*f2(x1+dx11/2., x2+dx12/2., t+dt/2.);
  
  dx31=dt*f1(x1+dx21/2., x2+dx22/2., t+dt/2.);
  dx32=dt*f2(x1+dx21/2., x2+dx22/2., t+dt/2.);
    
  dx41=dt*f1(x1+dx31, x2+dx32, t+dt);
  dx42=dt*f2(x1+dx31, x2+dx32, t+dt);
  
  x1 += dx11/6. + dx21/3. + dx31/3. + dx41/6.;
  x2 += dx12/6. + dx22/3. + dx32/3. + dx42/6.;
}



int main(void){

  double x1=0, x2=1, t=0;
  double dt=0.1;

  for(t=0; t<10; t+=dt){
    cout<<t<<" "<<x1<<" "<<fExacta(t)<<endl;
    UnPasoDeRungeKutta4_2(x1, x2, t, dt);
  }
    
  
  return 0;

}
