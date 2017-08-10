//Mi primer programa en C++


//Ctrl-k corta toda la linea a partir de ahi y Ctrl-y lo pega

#include <iostream> //print
#include <fstream> //save to file
#include <cmath> //math functions

using namespace std;

const double ERR = 1e-7;


double f(double x){
  return sin(x)/x;

}

double calcularCero(double a, double b){

//calcular el cero de f(x)
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
  cout<<"y el error es "<<m-M_PI<<endl; 


}

int main(void){

  //cout<<pow(i,2)<<endl;

  /*for(int i=0; i<10; cout<<i<<endl,i++ ){//mande a la pantalla i y luego un end line
 
    }
  int i, suma;
  for(suma=0,i=0;i<10;i++){
    suma+=2*i+1;
    

    }cout<<suma<<endl;

  for(i=0;i<10;i++){
    
    cout<<i<<endl;
    }


  
  double x;
  for(x=0.1; x<10; x+=0.1)
    cout<<x<<" "<<f(x)<<endl;
  



  //calcular el cero de f(x)
  double a,b,m;
  double fa,fb,fm;

  a=2; b=4;

  while(abs(b-a)>ERR){
    
    m=(a+b)/2.;
  if(f(a)*f(m)<0){
    b=m;
  }else{
    a=m;
  }
  
  }


  cout<<"El cero esta en x="<<m<<endl;
  cout<<"y el error es "<<m-M_PI<<endl;  */

  calcularCero(2, 4);
  
  return 0;

}
