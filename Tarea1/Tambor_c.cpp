#include <iostream>
#include <fstream> 
#include <cmath>
#include <sstream>

using namespace std;

const double alpha=0;
const double ERR = 1e-7;
const double r_a=5;

double f(double alpha, double x, double t){
  return cos(alpha*t-x*sin(t));
}


double IntegralPorSimpson(double alpha, double x, double a, double b, int N){

  double h=(b-a)/N, suma, t;
  int i;

  for(suma=0, i=0; i<=N; i++){

    t=a+i*h;
    if(i==0 || i==N){
      suma+=f(alpha, x, t);
    }else if(i%2==0){
      suma+=2*f(alpha, x, t);
    }else{
      suma+=4*f(alpha, x, t);
    }

  }
  return suma*h/3;
}


double Bessel(double alpha, double x){
  return 1/M_PI*IntegralPorSimpson(alpha, x, 0, M_PI, 100);
}


double f(double lambda){
  double r=r_a;
  double x=r*lambda;

  return Bessel(alpha, x);
}


double CeroBiseccion(double a, double b){

  double m;
  while(abs(b-a)>ERR){
    
  m=(a+b)/2.;
  if(f(a)*f(m)<0){
    b=m;
  }else{
    a=m;
  }
  
  }
  
  return m;

}



int main(void){

  double n_a,n_b,m;
  for(n_a=0.167, n_b=0.793; n_b<3.923; n_b+=0.626){

    double lambda=CeroBiseccion(n_a,n_b);
    cout<<"R(r=5)=0 en lambda="<<lambda<<endl;
    n_a=n_b;

    ofstream R_r;
    R_r.open("R_r"+std::to_string(lambda)+".txt");
    
    for(double r=0.01;r<=5.0;r+=0.01){
      R_r << r <<" "<< Bessel(alpha, r*lambda) << "\n";
    }
    R_r.close();
  }

  
  return 0;

}
