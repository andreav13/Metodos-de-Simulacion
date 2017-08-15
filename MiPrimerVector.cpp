#include <iostream>
#include <fstream> //save to file
#include <cmath> //math functions
#include "Vector.h"

using namespace std;


int main(void){

  vector3D a,b,c,d;

  /*a.cargue(1,2,3);
  b.cargue(-1,-2,-3);
  cout<<a.x()<<endl;
  a.show();
  c.operator=(a.operator+(b));
  c=a+b;
  c.show();
  b=a.operator*=(3);
  b=a*3; //b=a.operator*(3)
  b.show(); 
  b=3*a; //b=operator*(3,a)
  b.show();*/

  a.cargue(1,0,0);
  b.cargue(0,1,0);
  //c=b^a;           //producto cruz
  //c.show();

  c.cargue(0,0,1);
  d=a^b+c;  //no da (0,0,2) como esperado porque order of operations: MDAS luego operaciones logicas
  d=(a^b)+c;
  d.show();
  
  return 0;

}


