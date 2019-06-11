#include <iostream>
#include <cmath>
#include <stdio.h> 
#include "Random64.h"
using namespace std;

// const int Lx=256,Ly=256;	//Lx=200
const int Lx=256,Ly=256;	//Lx=200
const int Q=4;
const double p=0.25, p0=0.25;

double gauss(double x, double mu, double sigma){
  double g, aux1 = sigma*sqrt(2*M_PI), aux2=M_PI*(x-mu)*(x-mu)/(aux1*aux1);
  g = (1.0/aux1)*exp(-aux2);
  return g;
}

double tripleIntegral(double n[Lx][Ly][Q]){
  double ax[Lx], ay[Lx][Ly], answer=0.0, untercio=(1./3.);
  
  for(int i=0;i<Lx;i++)
    for(int j=0;j<Ly;j++){
      ay[i][j]=0;
      for(int k=0;k<Q;k++){
	ay[i][j]+=n[i][j][k];
      }
    }
  
  for(int i=0;i<Lx;i++){ 
    ax[i] = 0;
    for(int j=0;j<Ly;j++){ 
      if(j==0||j==Ly-1) 
	ax[i]+=ay[i][j]; 
      else if (j%2==0) 
	ax[i]+=2*ay[i][j]; 
      else
	ax[i]+=4*ay[i][j]; 
    }
    ax[i]=untercio*ax[i];
  }
  
  for(int i=0; i<Lx;i++){
    if(i==0||i==Lx-1) 
      answer+=ax[i]; 
    else if (i%2==0) 
      answer+=2*ax[i]; 
    else
      answer+=4*ax[i]; 
  } 
  answer*=untercio; 
  
  return answer; 
} 

class LatticeGas{
private:
  int V[2][Q];
  double n[Lx][Ly][Q], nnew[Lx][Ly][Q];
public:
  void Inicie(double mux,double muy,double sigmax,double sigmay);
  void Show(bool ImpNew);
  void Colisione(void);
  void Adveccione(void);
  double GetSigma2(void);
  double Integrese(void);
};
void LatticeGas::Inicie(double mux,double muy,double sigmax,double sigmay){
  //Definir los vectores
  V[0][0]=1;	V[0][1]=0;	V[0][2]=-1;	V[0][3]=0;
  V[1][0]=0;	V[1][1]=1;	V[1][2]=0;	V[1][3]=-1;
  int i, j, k;
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Q;k++){
	n[i][j][k]=0.25*gauss(i,mux,sigmax)*gauss(j,muy,sigmay);	//Asigna valores iniciales
	nnew[i][j][k]=0;	//Borra todos
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
void LatticeGas::Colisione(void){
  double aux=1-2*p-p0;
  for(int i=0;i<Lx;i++){
    for(int j=0;j<Ly;j++){
      for(int k=0;k<Q;k++)
	nnew[i][j][k]= p0*n[i][j][k]+p*n[i][j][(k+1)%Q]+p*n[i][j][(k+3)%Q]+aux*n[i][j][(k+2)%Q];
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
  long double N,r,miu,sigma2;int i,j,k;
  double z[Lx][Ly][Q];
    
    // Calculo miu
    for(i=0;i<Lx;i++)
      for(j=0;j<Ly;j++){
		r=sqrt(i*i+j*j);
		for(k=0;k<Q;k++){
		  z[i][j][k]=r*n[i][j][k];
		}
      }
  miu=tripleIntegral(z);
  
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++){
      r=sqrt(i*i+j*j);
      for(k=0;k<Q;k++){
	z[i][j][k]=(r-miu)*(r-miu)*n[i][j][k];
      }
    }
	sigma2 = tripleIntegral(z);

  return sigma2;
}


int main(void){
  LatticeGas Difusion;
  bool ImpNew=false;
  double mux=Lx/2.0,sigmax=8,muy=Ly/2.0,sigmay=8;
  
  Difusion.Inicie(mux,muy,sigmax,sigmay);
  // Difusion.Show(ImpNew);
  for(int t=0;t<350;t++){
    cout<<t<<" "<<Difusion.GetSigma2()<<endl;
    // Difusion.GetSigma2();
    Difusion.Colisione();
    Difusion.Adveccione();
    // Difusion.Show(ImpNew);
  }
  
  return 0;
}
