#include <iostream> //print
#include <fstream> //save to file
#include <cmath> //math functions

using namespace std;
const double ERR=1E-4;

double fExacta(double Omega0, double t){
  return sin(Omega0*t); 
}


double f1(double x1, double x2, double t){
  return x2;
}


double f2(double Omega0, double x1, double x2, double t){
  return -Omega0*Omega0*x1;
}


void UnPasoDeRungeKutta4_2(double Omega0, double &x1, double &x2, double &t, double dt){
  double dx11, dx21, dx31, dx41;
  double dx12, dx22, dx32, dx42;
  
  dx11=dt*f1(x1, x2, t);
  dx12=dt*f2(Omega0, x1, x2, t);
  
  dx21=dt*f1(x1+dx11/2., x2+dx12/2., t+dt/2.);
  dx22=dt*f2(Omega0, x1+dx11/2., x2+dx12/2., t+dt/2.);
  
  dx31=dt*f1(x1+dx21/2., x2+dx22/2., t+dt/2.);
  dx32=dt*f2(Omega0, x1+dx21/2., x2+dx22/2., t+dt/2.);
    
  dx41=dt*f1(x1+dx31, x2+dx32, t+dt);
  dx42=dt*f2(Omega0, x1+dx31, x2+dx32, t+dt);
  
  x1 += dx11/6. + dx21/3. + dx31/3. + dx41/6.;
  x2 += dx12/6. + dx22/3. + dx32/3. + dx42/6.;
}


double f(double Omega0){
  double x1=0, x2=1, t=0;
  double dt=0.01;
  for(t=0; t<1; t+=dt){
    UnPasoDeRungeKutta4_2(Omega0, x1, x2, t, dt);
  }
  return x1;
}


double calcularCero(double a, double b){
  double m;
  while(abs(b-a)>ERR){
    
  m=(a+b)/2.;
  if(f(a)*f(m)<0){
    b=m;
  }else{
    a=m;
  }
  }

  cout<<"El cero esta en x="<<m<<endl;

  }


int main(void){

  /*double Omega0;
  for(Omega0=0.1; Omega0<20; Omega0+=0.1){
    cout<<Omega0<<" "<<f(Omega0)<<endl;
    }*/
    
  calcularCero(2, 5);
  
  return 0;

}
