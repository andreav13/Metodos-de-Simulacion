#include <iostream> //print
#include <fstream> //save to file
#include <cmath> //math functions
using namespace std;


const double ERR = 1e-7;


double f(double alpha, double x, double t){
  return cos(alpha*t-x*sin(t));
}


double IntegralPorSimpson(double alpha, double x, double a, double b, int N){

  double h=(b-a)/N, suma, t;  //N tiene que ser par
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
  return 1/M_PI*IntegralPorSimpson(alpha, x, 0, M_PI, 50);
}



double calcularCero(double a, double b){

  double m, alpha=0;
  while(abs(b-a)>ERR){
    
  m=(a+b)/2.;
  if(Bessel(alpha,a)*Bessel(alpha,m)<0){
    b=m;
  }else{
    a=m;
  }
  
  }
  
  cout<<"El cero esta en t="<<m<<endl;

}



int main(void){

  calcularCero(0,5);
  calcularCero(5,7);
  calcularCero(7,10);
  
  return 0;

}
