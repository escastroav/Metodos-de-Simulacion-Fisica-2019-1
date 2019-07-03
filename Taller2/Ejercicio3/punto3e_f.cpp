#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=512;
const int Ly=64;

const int Q=9;

const double tau=1.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
const double Eta = (tau-0.5)/3.0;

const double Dt=1.0;

const int R = Ly/8;//El radio del círculo
const int R2 = R*R;
const int ixc=Lx/4;//Coordenadas del centro del círculo
const int iyc=Ly/2;
const int N =24;//El número de caras del polígono regular en que se divide la circunferencia
const double A = 2.0*M_PI*R/N;//"Area" de una cara del polígono
const double PI2= 2*M_PI;
const double PI2oN=PI2/N;

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
  double sigmaxx(int ix,int iy,bool UseNew);
  double sigmaxy(int ix,int iy,bool UseNew);
  double sigmayy(int ix,int iy,bool UseNew);
  double Fx(double x, double y, double Ax, double Ay);
  double Fy(double x, double y, double Ax, double Ay);
  double Ax(double Theta);//Componentes del area del lado
  double Ay(double Theta);
  double Px(double Theta);//Punto donde quiero medir la fuerza
  double Py(double Theta);
  double feq(double rho0,double Ux0,double Uy0,int i);
  void Colisione(void);
  void Adveccione(void);
  void Inicie(double rho0,double Ux0,double Uy0);
  void ImponerCampos(double Vventilador, double omega);
  void Imprimase(const char * NombreArchivo,double Vventilador);
  void TotalF(double Vventilador, double omega);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0]=4/9.0; 
  w[1]=w[2]=w[3]=w[4]=1/9.0;
  w[5]=w[6]=w[7]=w[8]=1/36.0;
  //Cargar los vectores
  V[0][0]=0;
  V[1][0]=0;

  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1; V[0][4]=0;
  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;

  V[0][5]=1;  V[0][6]=-1; V[0][7]=-1; V[0][8]=1;
  V[1][5]=1;  V[1][6]=1;  V[1][7]=-1; V[1][8]=-1;
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

//calculo del tensor de esfuerzos

double LatticeBoltzmann::sigmaxx(int ix,int iy,bool UseNew){
  double p = rho(ix, iy, UseNew)/3.0;
  double dUxdx;int i;
  for(dUxdx=0,i=0;i<Q;i++)
    dUxdx+=w[i]*V[0][i]*Jx((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly, UseNew)/rho((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly, UseNew);
  return -p+6*Eta*dUxdx/Dt;
}

double LatticeBoltzmann::sigmaxy(int ix,int iy,bool UseNew){
  double dUxdy, dUydx;int i;
  for(dUxdy=0,dUydx=0,i=0;i<Q;i++){
    dUxdy+=w[i]*V[1][i]*Jx((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly, UseNew)/rho((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly, UseNew);
    dUydx+=w[i]*V[0][i]*Jy((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly, UseNew)/rho((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly, UseNew);
  }
}

double LatticeBoltzmann::sigmayy(int ix,int iy,bool UseNew){
  double p = rho(ix, iy, UseNew)/3.0;
  double dUydy;int i;
  for(dUydy=0,i=0;i<Q;i++)
    dUydy+=w[i]*V[1][i]*Jy((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly, UseNew)/rho((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly, UseNew);
  return -p+6*Eta*dUydy/Dt;
}

double LatticeBoltzmann::Fx(double x, double y, double Ax, double Ay){
  int ix = (int) x; int iy = (int) y;
  double u = (x-ix), v =(y-iy);
  double sigmaxx0, sigmaxy0;
  bool UN=false;
  //Interpolación del tensor de esfuerzos, sigma=F/A
  sigmaxx0=sigmaxx(ix,iy,UN)*(1-u)*(1-v)+sigmaxx(ix+1,iy,UN)*u*(1-v)+sigmaxx(ix,iy+1,UN)*(1-u)*v+sigmaxx(ix+1,iy+1,UN)*u*v;
  sigmaxy0=sigmaxy(ix,iy,UN)*(1-u)*(1-v)+sigmaxy(ix+1,iy,UN)*u*(1-v)+sigmaxy(ix,iy+1,UN)*(1-u)*v+sigmaxy(ix+1,iy+1,UN)*u*v;
  return sigmaxx0*Ax+sigmaxy0*Ay;
}

double LatticeBoltzmann::Fy(double x, double y, double Ax, double Ay){
  int ix = (int) x; int iy = (int) y;
  double u = (x-ix), v =(y-iy);
  double sigmayy0, sigmaxy0;
  bool UN=false;
  //Interpolación del tensor de esfuerzos, sigma=F/A
  sigmayy0=sigmayy(ix,iy,UN)*(1-u)*(1-v)+sigmayy(ix+1,iy,UN)*u*(1-v)+sigmayy(ix,iy+1,UN)*(1-u)*v+sigmayy(ix+1,iy+1,UN)*u*v;
  sigmaxy0=sigmaxy(ix,iy,UN)*(1-u)*(1-v)+sigmaxy(ix+1,iy,UN)*u*(1-v)+sigmaxy(ix,iy+1,UN)*(1-u)*v+sigmaxy(ix+1,iy+1,UN)*u*v;
  return sigmayy0*Ay+sigmaxy0*Ax;
}

double LatticeBoltzmann::Ax(double Theta){
  return A*cos(Theta);
}

double LatticeBoltzmann::Ay(double Theta){
  return A*sin(Theta);
}

double LatticeBoltzmann::Py(double Theta){
  return R*sin(Theta)+1.0*iyc;
}

double LatticeBoltzmann::Px(double Theta){
  return R*cos(Theta)+1.0*ixc;
}

double LatticeBoltzmann::feq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}
void LatticeBoltzmann::Colisione(void){
  int ix,iy,i; double rho0,Ux0,Uy0;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macroscópicas
      rho0=rho(ix,iy,false);  Ux0=Jx(ix,iy,false)/rho0;  Uy0=Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++)
	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Ux0,Uy0,i);
    }
}
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}

void LatticeBoltzmann::Inicie(double rho0,double Ux0,double Uy0){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[ix][iy][i]=feq(rho0,Ux0,Uy0,i);
}

void LatticeBoltzmann::ImponerCampos(double Vventilador, double omega){
  int i,ix,iy; double rho0, Ux0, Uy0;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      //ventilador
      if(ix==0) 
	for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,Vventilador,0,i);
      //Obstáculo cilíndrico
      else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R2){
	Ux0=-omega*(iy-iyc); Uy0=omega*(ix-ixc);//poner a rotar al cilindro con velocidad angular omega
	for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,Ux0,Uy0,i);}
    }
}

void LatticeBoltzmann::Imprimase(const char * NombreArchivo,double Vventilador){
  ofstream MiArchivo(NombreArchivo); double rho0,Ux0,Uy0;
  for(int ix=0;ix<Lx;ix+=4){
    for(int iy=0;iy<Ly;iy+=4){
      rho0=rho(ix,iy,true);  Ux0=Jx(ix,iy,true)/rho0;  Uy0=Jy(ix,iy,true)/rho0;
      MiArchivo<<ix<<" "<<iy<<" "<<4*(Ux0-Vventilador)/Vventilador<<" "<<4*Uy0/Vventilador<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}

void LatticeBoltzmann::TotalF(double Vventilador, double omega){
  double FxT=0, FyT=0;
  double Axk,Ayk, Pxk,Pyk, Theta, rho0;
  for(int k=0;k<N;k++){
    Theta=k*PI2oN;
    Axk=Ax(Theta); Ayk=Ay(Theta); Pxk=Px(Theta); Pyk=Py(Theta);
    FxT+=Fx(Pxk,Pyk,Axk,Ayk); FyT+=Fy(Pxk,Pyk,Axk,Ayk);
  }
  rho0=rho(200,32,false);
  std::cout<<omega<<"\t"<<-2*omega*Vventilador*rho0*R*R<<"\t"<<FxT<<"\t"<<FyT<<"\n";
}

int main(void){
  LatticeBoltzmann Aire;
  int t,tmax=1000;
  double RHOinicial=1.0, Vventilador=0.1;

  double omegaIni=PI2/1000, omegaMax=omegaIni+1e-3, domega=(omegaMax-omegaIni)/20.0;
  Aire.Inicie(RHOinicial,Vventilador,0);
  for(double omega=omegaIni;omega<=omegaMax;omega+=domega){
    for(t=0;t<tmax;t++){
      Aire.Colisione();
      Aire.ImponerCampos(Vventilador, omega);
      Aire.Adveccione();
    }
    
    Aire.TotalF(Vventilador, omega);//Para magnus la que importa es Fy, F=omegaxVvent y omega va en z y Vvent en x
  }
  return 0;
}
