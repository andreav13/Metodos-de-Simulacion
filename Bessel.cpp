#include <iostream> //print
#include <fstream> //save to file
#include <cmath> //math functions

using namespace std;


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



int main(void){
  double x, alpha=0;
  for(x=0;x<20;x+=0.1){
    cout<<x<<" "<<Bessel(alpha, x)<<endl;
  }

  

  //cout<<"Integral="<<IntegralPorSimpson(0,M_PI,10)<<endl;
  
  return 0;

}
