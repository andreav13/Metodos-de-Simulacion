#include <iostream> //print
#include <fstream> //save to file
#include <cmath> //math functions

using namespace std;

double f(double x){
  return sin(x);

}

double IntegralPorSimpson(double a, double b, int N){

  double h=(b-a)/N, suma, x;  //N tiene que ser par
  int i;

  for(suma=0, i=0; i<=N; i++){

    x=a+i*h;
    if(i==0 || i==N){
      suma+=f(x);
    }else if(i%2==0){
      suma+=2*f(x);
    }else{
      suma+=4*f(x);
    }

  }
  return suma*h/3;

}






int main(void){


  int N; double h, a=0, b=M_PI;
  for(N=2;N<1025;N*=2){
    h=(b-a)/N;
    cout<<h<<" "<<IntegralPorSimpson(a,b,N)-2<<endl;
    


  }

  

  //cout<<"Integral="<<IntegralPorSimpson(0,M_PI,10)<<endl;
  
  return 0;

}
