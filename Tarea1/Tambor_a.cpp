#include <iostream>
#include <fstream>
#include <cmath> 

using namespace std;

const double lambda=1;


double f1(double R, double x2, double r){
  return x2;
}


double f2(double R, double x2, double r){
  return -1/r*x2-lambda*lambda*R;
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

  double R=1, x2=0, r;
  double dr=0.01;

  for(r=0.01; r<10.0; r+=dr){
    cout<<r<<" "<<R<<endl;
    UnPasoDeRungeKutta4_2(R, x2, r, dr);
  }
    
  
  return 0;

}
