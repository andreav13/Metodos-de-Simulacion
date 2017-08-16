#include <iostream>
#include <fstream> 
#include <cmath> 

using namespace std;

const double B=0.3;
const double Y=0.05;


double f_s(double s, double i, double t){
  return -B*s*i;
}


double f_i(double s, double i, double t){
  return B*s*i-Y*i;
}

void UnPasoDeRungeKutta4_2(double &s, double &i, double &r, double &t, double dt){
  double ds1, ds2, ds3, ds4;
  double di1, di2, di3, di4;
  
  ds1=dt*f_s(s, i, t);
  di1=dt*f_i(s, i, t);
  
  ds2=dt*f_s(s+ds1/2., i+di1/2., t+dt/2.);
  di2=dt*f_i(s+ds1/2., i+di1/2., t+dt/2.);
  
  ds3=dt*f_s(s+ds2/2., i+di2/2., t+dt/2.);
  di3=dt*f_i(s+ds2/2., i+di2/2., t+dt/2.);
    
  ds4=dt*f_s(s+ds3, i+di3, t+dt);
  di4=dt*f_i(s+ds3, i+di3, t+dt);
  
  s += ds1/6. + ds2/3. + ds3/3. + ds4/6.;
  i += di1/6. + di2/3. + di3/3. + di4/6.;
  r = 1-s-i;
}




int main(void){

  double s=0.999, i=0.001, r=1-s-i, t=0;
  double dt=.5;

  for(t=0; t<70; t+=dt){
    cout<<t<<" "<<s<<" "<<i<<" "<<r<<endl;
    UnPasoDeRungeKutta4_2(s, i, r, t, dt);
  }
  

  return 0;
}
