#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;

// const int Lx=256,Ly=256;	//Lx=200
const int Lx=256,Ly=256;	//Lx=200
const int Q=4;
const double p=0.25, p0=0.25;

class LatticeGas{
private:
  int V[2][Q];
  int n[Lx][Ly][Q], nnew[Lx][Ly][Q];
public:
  void Inicie(int Bolitas,double mux,double muy,double sigmax,double sigmay,Crandom & ran64);
  void Show(bool ImpNew);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double GetSigma2(void);
};
void LatticeGas::Inicie(int Bolitas,double mux,double muy,double sigmax,double sigmay,Crandom & ran64){
  //Definir los vectores
  V[0][0]=1;	V[0][1]=0;	V[0][2]=-1;	V[0][3]=0;
  V[1][0]=0;	V[1][1]=1;	V[1][2]=0;	V[1][3]=-1;
  int i, j, k, b;
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Q;k++){
	n[i][j][k]=0;	//Borra todos
	nnew[i][j][k]=0;	//Borra todos
      }  
  while(Bolitas>0){	//Coloca Bolitas
    i=(int) ran64.gauss(mux,sigmax);
    j=(int) ran64.gauss(muy,sigmay);
    if(j<0) j=0;
    if(i<0) i=0;
    if(i>Lx-1) i=Lx-1;
    if(j>Ly-1) j=Ly-1;
    k=(int) Q*ran64.r();
    if(n[i][j][k]==0){
      n[i][j][k]=1;
      Bolitas--;
    }
  }
}
void LatticeGas::Show(bool ImpNew){
  int i, j, k;
  if(ImpNew){
    for(k=0;k<Q;k++){
      for(j=0;j<Ly;j++){
	for(i=0;i<Lx;i++)
	  cout<<nnew[i][j][k];
	cout<<endl;
      }
      cout<<endl;
    }
  } 
  else{
    for(k=0;k<Q;k++){
      for(j=0;j<Ly;j++){
	for(i=0;i<Lx;i++)
	  cout<<n[i][j][k];
	cout<<endl;
      }
      cout<<endl;
    }	
  }
  cout<<endl;
}
void LatticeGas::Colisione(Crandom & ran64){
  double aux;
  for(int i=0;i<Lx;i++){
    for(int j=0;j<Ly;j++){
      aux=ran64.r();
      if(aux<=p0){
	for(int k=0;k<Q;k++){
	  for(int l=0;l<2;l++){
	    nnew[i][j][k]=n[i][j][k];  
	  } 
	}
      }
      else if(aux>=p0 && aux<p+p0){
	for(int k=0;k<Q;k++)
	  nnew[i][j][(k+1)%Q]=n[i][j][k];
      }
      else if(aux>=p+p0 && aux<2*p+p0){
	for(int k=0;k<Q;k++)
	  nnew[i][j][(k-1)%Q]=n[i][j][k];
      }
      else{
	for(int k=0;k<Q;k++)
	  nnew[i][j][(k+2)%Q]=n[i][j][k];
      }
    } 
  }
}
void LatticeGas::Adveccione(void){
  for(int i=0;i<Lx;i++){		//Para cada celda
    for(int j=0;j<Ly;j++){
      for(int k=0;k<Q;k++){	//Para cada direccion
	  n[(i+Lx+V[0][k])%Lx][(j+Ly+V[1][k])%Ly][k]=nnew[i][j][k];		//La celda aledanha en direccion V[k][l]
      }	
    }
  }
}
double LatticeGas::GetSigma2(void){
  double N,r,rprom,sigma2;int i,j,k;
  // Calculo N
  for(N=0,i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Q;k++)
	  N+=n[i][j][k];
  // Calculo rprom
  for(rprom=0,i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Q;k++){
	  rprom+=sqrt(i*i+j*j)*n[i][j][k];
	}
  rprom/=N;
  // Calculo sigma2
  for(sigma2=0,i=0;i<Lx;i++)
    for(j=0;j<Ly;j++){
		r=sqrt(i*i+j*j);
		for(k=0;k<Q;k++){
		  sigma2+=r*n[i][j][k]*(r*n[i][j][k]-2*rprom);
	  }
	}
	sigma2+=N*rprom*rprom;
  sigma2/=(N-1);
  return sigma2;
}

int main(void){
  LatticeGas Difusion;
  Crandom ran64(1);
  bool ImpNew=false;
  int Bolitas=1200;
  double mux=Lx/2.0,sigmax=16,muy=Ly/2.0,sigmay=16;
  
  Difusion.Inicie(Bolitas,mux,muy,sigmax,sigmay,ran64);
  // Difusion.Show(ImpNew);
  // cout<<endl;
  for(int t=0;t<350;t++){
    cout<<t<<" "<<Difusion.GetSigma2()<<endl;
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
    // Difusion.Show(ImpNew);
  }
  
  return 0;
}
