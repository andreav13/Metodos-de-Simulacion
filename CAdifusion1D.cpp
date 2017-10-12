#include <iostream>
#include <cmath> 
#include "Random64.h"

using namespace std;

const int Lx=200;
const double p=0.5;

class LatticeGas{
private:
  int V[2]; //V[i] i=0:derecha, i=1:izquierda
  int n[Lx][2], nnew[Lx][2]; //n[ix][i]
public:
  LatticeGas(void); //constructor
  void Inicie(Crandom & ran64);
  void Show(int copia);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  int DondeEstaLaBolita(void);
};

LatticeGas::LatticeGas(void){
  V[0]=1; V[1]=-1;
}

void LatticeGas::Inicie(Crandom & ran64){
  for(int ix=0;ix<Lx;ix++)
    for(int i=0;i<2;i++){
      n[ix][i]=0; nnew[ix][i]=0;
    }
  if(ran64.r()<0.5)
       n[Lx/2][0]=1;
  else
    n[Lx/2][1]=1;
}

void LatticeGas::Show(int copia){
  for(int i=0;i<2;i++){
    for(int ix=0;ix<Lx;ix++){
      if(copia==0)
	cout<<n[ix][i];
      else
	cout<<nnew[ix][i];
    }
    cout<<endl;
  }
}

void LatticeGas::Colisione(Crandom & ran64){ //de n a nnew
  int ix;
  for(ix=0;ix<Lx;ix++){
    if(ran64.r()<p){
      nnew[ix][0]=n[ix][0]; nnew[ix][1]=n[ix][1]; 
    }
    else{
      nnew[ix][0]=n[ix][1]; nnew[ix][1]=n[ix][0];
    }
  }
}
   
void LatticeGas::Adveccione(void){ //de nnew a n
  int ix,i;
  for(ix=0;ix<Lx;ix++)
    for(i=0;i<2;i++)
      n[(ix+V[i]+Lx)%Lx][i]=nnew[ix][i]; //para que el ultimo lleve al primero
                                         //(condiciones de frontera periodicas)
}

int LatticeGas::DondeEstaLaBolita(void){
  int ix=0;
  while(n[ix][0]+n[ix][1]==0)
    ix++;
  return ix;
}

//-----------------------Funciones Globales-----------------------

const int N=1000;

double sigma2(LatticeGas * Difusion){
  double Xprom,sigma2; int c;
  //calculo xprom
  for(Xprom=0,c=0;c<N;c++)
    Xprom+=Difusion[c].DondeEstaLaBolita();
  Xprom/=N;
  //calculo sigma2
  for(sigma2=0,c=0;c<N;c++)
    sigma2+=pow(Difusion[c].DondeEstaLaBolita()-Xprom, 2);
  sigma2/=(N-1);
  
  return sigma2;
}


int main(void){

  Crandom ran64(1012);
  LatticeGas Difusion[N];
  int c,t,tmax=400;

  //Inicie
  for(c=0;c<N;c++) Difusion[c].Inicie(ran64);

  //Corra
  for(t=0;t<tmax;t++){
    cout<<t<<" "<<sigma2(Difusion)<<endl;
    for(c=0;c<N;c++) Difusion[c].Colisione(ran64);
    for(c=0;c<N;c++) Difusion[c].Adveccione();
  }


  /*Difusion.Show(0); //show n
  Difusion.Colisione(ran64);
  Difusion.Show(1); //show nnew
  Difusion.Adveccione();
  Difusion.Show(0);*/
  
  return 0;
}
