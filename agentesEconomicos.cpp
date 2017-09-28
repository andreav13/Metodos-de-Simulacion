#include <iostream>
#include <cmath>
#include "Random64.h"

using namespace std;

const int N=10;

class Agente{
private:
  double dinero;
public:
  void Inicie(double dinero0){dinero=dinero0;};
  double GetDinero(void){return dinero;};

  
  friend class Mercado;

};

  
class Mercado{
private:

public:
  void HagaTransaccionEntre(Agente & Agente1, Agente & Agente2, Crandom & ran64);

};


void Mercado::HagaTransaccionEntre(Agente & Agente1, Agente & Agente2, Crandom & ran64){
  double GananciaAgente1, GananciaAgente2, Resultado;
  Resultado=2000*ran64.r();
  GananciaAgente1 = Resultado-1000;
  Agente1.dinero+=GananciaAgente1;
  GananciaAgente2 = 2000-Resultado-1000;
  Agente2.dinero+=GananciaAgente2;
}


  

int main(void){

  Mercado Corabastos;
  Agente Paisano[N];  
  Crandom ran64(1);
  int i,j,t,Ntransacciones=1;

  //Inicie todos los agentes con 10000 pesos
  for(i=0;i<N;i++)
    Paisano[i].Inicie(10000);
  
  for(t=0;t<Ntransacciones;t++){

    //Escojo dos agentes al azar
    i=(int) N*ran64.r();
    j=(int) (N-1)*ran64.r(); if(j>i) j++;
    
    //Los pongo a interactuar
    Corabastos.HagaTransaccionEntre(Paisano[i], Paisano[j], ran64);
  }
  
  for(i=0;i<N;i++)
    cout<<Paisano[i].GetDinero()<<endl;
  
  
  return 0;
}
