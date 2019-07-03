#include<iostream>
#include<cmath>
#include "Random64.h"

const int Lx=256;//Tamaño del lattice
const int Ly=256;
const double p0=0.25;//Probabilidad de que se quede quieto.
const double p=0.25;//Probabilidad girar a la derecha
const int Q = 4;//El número de velocidades
const int D = 2;//Número de dimensiones

double gaussiana(double x, double mu, double sigma){
  double aux1=pow(2.0*M_PI*sigma*sigma,-0.5), aux2=M_PI*aux1*aux1;
  return aux1*exp(-(x-mu)*(x-mu)*aux2);
}

double tripleIntegral(double n[Lx][Ly][Q]){
  double ax[Lx], ay[Lx][Ly], answer=0.0, untercio=1./3.0;
  
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
  int V[D][Q];//Velocidades
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q];//La matríz de celdas vectores, así como un nnew sobre el cuál escribir, f[i][j][0]=east, f[1]=north, f[2]=west, f[3]=south.
public:
  void Inicie(double mux, double sigmax, double muy, double sigmay);
  void Show(bool ImprimirNew);
  void Colisione();
  void Adveccione(void);
  double GetSigma2(void);
};

void LatticeGas::Inicie(double mux, double sigmax,  double muy, double sigmay){
  //Definimos los vectores
  V[0][0]=1;V[1][0]=0;
  V[0][1]=0;V[1][1]=1;
  V[0][2]=-1;V[1][2]=0;
  V[0][3]=0;V[1][3]=-1;
  //Llenamos los contenidos de las celdas
  int i,j,k;//Contadores de matríz, y direccion
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Q;k++)
	f[i][j][k]=fnew[i][j][k]=0.25*gaussiana(i,mux,sigmax)*gaussiana(j,muy,sigmay);//Inicio como una bigaussiana
}

void LatticeGas::Show(bool ImprimirNew){
  for(int k=0;k<Q;k++){
    for(int i=0;i<Lx;i++)
      for(int j=0;j<Ly;j++)
	if(ImprimirNew) std::cout<<fnew[i][j][k]; else std::cout<<f[i][j][k];
    std::cout<<std::endl;
  }
}

void LatticeGas::Colisione(){
  for(int i=0;i<Lx;i++)
    for(int j=0;j<Ly;j++)
      for(int k=0;k<Q;k++)
	fnew[i][j][k]=p0*f[i][j][k]+p*(f[i][j][(k+1)%Q]+f[i][j][(k+3)%Q])+(1-2*p-p0)*f[i][j][(k+2)%Q];
      /*      //0->este, 1->norte, 2->oeste, 3->sur
      f[i][j][0]=p0*f[i][j][0]+p*(f[i][j][1]+f[i][j][3])+(1-2*p-p0)*f[i][j][2];
      f[i][j][1]=p0*f[i][j][1]+p*(f[i][j][0]+f[i][j][2])+(1-2*p-p0)*f[i][j][3];
      f[i][j][2]=p0*f[i][j][2]+p*(f[i][j][1]+f[i][j][3])+(1-2*p-p0)*f[i][j][0];
      f[i][j][3]=p0*f[i][j][3]+p*(f[i][j][0]+f[i][j][2])+(1-2*p-p0)*f[i][j][1];*/
}

void LatticeGas::Adveccione(){//Muevo cada resultado de celda en la siguiente
  for(int i=0;i<Lx;i++)
    for(int j=0;j<Ly;j++)
      for(int k=0;k<Q;k++)
	  f[(i+Lx+V[0][k])%Lx][(j+Ly+V[1][k])%Ly][k]=fnew[i][j][k];//Condiciones periódicas a derecha e izquierda
}

double LatticeGas::GetSigma2(void){
  double xprom, yprom, sigma2x, sigma2y, sigma2; int i,j,k;
  double r, r2;
  double val[Lx][Ly][Q];
  //xprom
  for(xprom=0,i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Q;k++)
	val[i][j][k]=i*f[i][j][k];
  xprom=tripleIntegral(val);
  //sigma2x
  for(sigma2x=0,i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Q;k++)
	val[i][j][k]=i*i*f[i][j][k];
  sigma2x=tripleIntegral(val)-xprom*xprom;
  /*   //yprom
  for(yprom=0,i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Q;k++)
       val[i][j][k]=j*f[i][j][k];
       yprom=tripleIntegral(val);
  //sigma2y
  for(sigma2y=0,i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Q;k++)
	val[i][j][k]=j*j*f[i][j][k];
	sigma2y=tripleIntegral(val)-yprom*yprom;

  //rprom
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++){
      r=sqrt(i*i+j*j);
      for(k=0;k<Q;k++)
	val[i][j][k]=r*f[i][j][k];
  rprom=tripleIntegral(val);
  //sigmar
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++){
    r2=i*i+j*j;
      for(k=0;k<Q;k++)
	val[i][j][k]=r2*f[i][j][k];
    }
    sigma2 = tripleIntegral(val)-rprom*rprom;*/
  return sigma2x;
}

int main(void){
  LatticeGas Difusion;
  double mux=Lx/2.0, sigmax=16.0, muy=Ly/2.0, sigmay=16.0;
  int t, tmax=350;

  Difusion.Inicie(mux,sigmax,muy,sigmay);
  //  Difusion.Show(false);
  for(t=0;t<tmax;t++){
    std::cout<<t<<" "<<Difusion.GetSigma2()<<"\n";
    Difusion.Colisione();
    Difusion.Adveccione();
    }
  
  return 0;
}
