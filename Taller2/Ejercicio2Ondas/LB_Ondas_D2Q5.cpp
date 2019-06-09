#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=400;
const int Ly=200;

const int Q=5;
const double W0=1.0/3;

const double TanTheta = tan(M_PI/3.0);
//const double C=0.5; // C<0.707 celdas/click
//const double TresC2=3*C*C;
//const double AUX0=1-TresC2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  double feq(double rho0,double Jx0,double Jy0,int ix, int iy, int i);
  void Colisione(void);
  void Adveccione(void);
  void Inicie(double rho0,double Jx0,double Jy0);
  void ImponerCampos(int t);
  void Imprimase(const char * NombreArchivo);
  double Ccelda(int ix, int iy);
  double TresC2(int ix, int iy)
  {return 3*Ccelda(ix,iy)*Ccelda(ix,iy);};
  double AUX0(int ix, int iy)
  {return 1-TresC2(ix,iy)*(1-W0);};
};

LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1-W0)/4.0;
  //Cargar los vectores
  V[0][0]=0;
  V[1][0]=0;

  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1; V[0][4]=0;
  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]; else suma+=f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[0][i]; else suma+=f[ix][iy][i]*V[0][i];
  return suma;
}
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[1][i]; else suma+=f[ix][iy][i]*V[1][i];
  return suma;
}
double LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,int ix, int iy, int i){
  if(i==0)
    return rho0*AUX0(ix,iy);
  else
    return w[i]*(TresC2(ix,iy)*rho0+3*(V[0][i]*Jx0+V[1][i]*Jy0));
}
void LatticeBoltzmann::Colisione(void){
  int ix,iy,i; double rho0,Jx0,Jy0;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macroscópicas
      rho0=rho(ix,iy,false);  Jx0=Jx(ix,iy,false);  Jy0=Jy(ix,iy,false);
      for(i=0;i<Q;i++)
	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Jx0,Jy0,ix,iy,i);
    }
}
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}
void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[ix][iy][i]=feq(rho0,Jx0,Jy0,ix,iy,i);
}
void LatticeBoltzmann::ImponerCampos(int t){
  int i,ix=0,iy; double lambda,omega,rho0,Jx0,Jy0;
  lambda=10; //omega=2*M_PI/lambda;
  /*rho0=10*sin(omega*t)*/; Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){
	omega=2*M_PI*Ccelda(ix,iy)/lambda;
	rho0=10*sin(omega*t);
	fnew[0][iy][i]=feq(rho0,Jx0,Jy0,ix,iy,i);}
}
void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  ofstream MiArchivo(NombreArchivo); double rho0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);
      MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}
double LatticeBoltzmann::Ccelda(int ix, int iy)
{
  return (0.25-0.5)*0.5*tanh(ix-(-TanTheta*iy+100*(1+TanTheta)))+(0.5+0.25)*0.5;
  //return 0.5;
}
int main(void){
  LatticeBoltzmann Ondas;
  int t,tmax=400;
  
  Ondas.Inicie(0,0,0);
  for(t=0;t<tmax;t++){
    Ondas.Colisione();
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
  }
  Ondas.Imprimase("Ondas.dat");

  return 0;
}
