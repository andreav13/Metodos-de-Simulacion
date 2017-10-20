#include <iostream>
#include <cmath> 
#include "Random64.h"

using namespace std;

const int Lx=200;
const double p=0.5;

class LatticeGas{
private:
  int V[2]; //V[i] i=0:derecha, i=1:izquierda
  double n[Lx][2], nnew[Lx][2]; //n[ix][i]
public:
  LatticeGas(void); //constructor
  void Inicie(double mu, double sigma);
  void ShowRho();
  void Colisione();
  void Adveccione(void);
  double GetRho(int ix){return n[ix][0]+n[ix][1];};
  double GetSigma2(void);
};

LatticeGas::LatticeGas(void){
  V[0]=1; V[1]=-1;
}

void LatticeGas::Inicie(double mu, double sigma){
  for(int ix=0;ix<Lx;ix++)
    n[ix][0]=n[ix][1]=1./(2*sigma*sqrt(2*M_PI))*exp(-0.5*pow((ix-mu)/sigma,2));
}

void LatticeGas::ShowRho(){
  for(int ix=0;ix<Lx;ix++)
    cout<<ix<<" "<<GetRho(ix)<<endl;
}

void LatticeGas::Colisione(){ //de n a nnew
  int ix;
  for(ix=0;ix<Lx;ix++){
      nnew[ix][0]=p*n[ix][0]+(1-p)*n[ix][1];
      nnew[ix][1]=p*n[ix][1]+(1-p)*n[ix][0];
  }
}
   
void LatticeGas::Adveccione(void){ //de nnew a n
  int ix,i;
  for(ix=0;ix<Lx;ix++)
    for(i=0;i<2;i++)
      n[(ix+V[i]+Lx)%Lx][i]=nnew[ix][i]; //para que el ultimo lleve al primero
                                         //(condiciones de frontera periodicas)
}

double LatticeGas::GetSigma2(void){
  double N,Xprom,sigma2; int ix;
  //calculo N
  for(N=0,ix=0;ix<Lx;ix++)
    N+=GetRho(ix);
  //calculo xprom
  for(Xprom=0,ix=0;ix<N;ix++)
    Xprom+=ix*GetRho(ix);
  Xprom/=N;
  //calculo sigma2
  for(sigma2=0,ix=0;ix<Lx;ix++)
    sigma2+=pow(ix-Xprom, 2)*GetRho(ix);
  sigma2/=N;
  
  return sigma2;
}


//-----------------------Funciones Globales-----------------------


int main(void){

  LatticeGas Difusion;
  int t,tmax=100;

  //Inicie
  Difusion.Inicie(Lx/2,Lx/48);

  //Corra
  for(t=0;t<tmax;t++){
    //cout<<t<<" "<<Difusion.GetSigma2()<<endl;
    Difusion.Colisione();
    Difusion.Adveccione();
  }

  //Mostrar resultado
  Difusion.ShowRho();

  
  return 0;
}
