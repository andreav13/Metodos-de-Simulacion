#include <iostream> //print
#include <fstream> //save to file
#include <cmath> //math functions

using namespace std;

const double alpha=0;

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
  double r=1;
  double x=r*lambda;

  return Bessel(alpha, x);
  
}



int main(void){

  double lambda;
  
  for(lambda=0.01;lambda<=80;lambda+=0.01){
    cout<<lambda<<" "<<f(lambda)<<endl;
  }

  
  return 0;

}
