#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx = 100;
const int Ly = 100;

const int Q = 5;
const double W0 = 1.0/3;

const double C = 0.5; //C<1/sqrt(2) o feq se hace 0 en i=0
const double TresC2 = 3*C*C;
const double AUX0=1-TresC2*(1-W0);

const double tau = 0.5;
const double Utau = 1.0/tau;
const double UmUtau = 1-Utau;

class LatticeBoltzmann
{
private:
  double w[Q];
  int V[2][Q]; //V[0][i]=V_ix ; V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q];
  public:
  LatticeBoltzmann(void);
  //Calcular cantidades macroscopicas
  double rho(int ix, int iy, bool useNew);
  double Jx(int ix, int iy, bool useNew);
  double Jy(int ix, int iy, bool useNew);
  double feq(double rho0, double Jx0, double Jy0, int i);
  //Colision y Advencion
  void Colisione(void);
  void Adveccione(void);
  //Inicializar
  void Inicie(double rho0, double Jx0, double Jy0);
  //Condicion De Contorno
  void ImponerCampos(int t);
  //Imprimir
  void Imprimase(const char * NombreArchivo);
};

LatticeBoltzmann::LatticeBoltzmann(void)
{
  //cargar los pesos
  w[0] = W0; w[1]=w[2]=w[3]=w[4]=(1-W0)/4.0;
  //Cargar los vectores
  V[0][0]=0;
  V[1][0]=0;

  V[0][1]=1;
  V[1][1]=0;

  V[0][2]=0;
  V[1][2]=1;

  V[0][3]=-1;
  V[1][3]=0;

  V[0][4]=0;
  V[1][4]=-1;
}
//Calcular cantidades macroscopicas
double LatticeBoltzmann::rho(int ix, int iy, bool useNew)
{
  double suma;
  for(int i=0,suma = 0;i<Q;i++)
    {
      if(useNew)
	suma+=fnew[ix][iy][i];
      else
	suma+=f[ix][iy][i];
    }

  return suma;
}
double LatticeBoltzmann::Jx(int ix, int iy, bool useNew)
{
  double suma;
  for(int i=0,suma = 0;i<Q;i++)
    {
      if(useNew)
	suma+=fnew[ix][iy][i]*V[0][i];
      else
	suma+=f[ix][iy][i]*V[0][i];
    }

  return suma;
}
double LatticeBoltzmann::Jy(int ix, int iy, bool useNew)
{
  double suma;
  for(int i=0,suma = 0;i<Q;i++)
    {
      if(useNew)
	suma+=fnew[ix][iy][i]*V[1][i];
      else
	suma+=f[ix][iy][i]*V[1][i];
    }

  return suma;
}

double LatticeBoltzmann::feq(double rho0, double Jx0, double Jy0, int i)
{
  if(i==0)
    return rho0*AUX0;
  else
    return w[i]*(TresC2*rho0+3*(V[0][i]*Jx0+V[1][i]*Jy0));
}
void LatticeBoltzmann::Colisione(void)
{
  double rho0=0, Jx0=0, Jy0=0;
  //Para cada celda
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      {
	//Calcular cantidades macroscopicas
	rho0=rho(ix,iy,false); Jx0=Jx(ix, iy, false); Jy0=Jy(ix,iy,false);
	//Colision
	for(int i=0; i<Q; i++)
	  fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0, Jx0, Jy0, i);
      }
}
void LatticeBoltzmann::Adveccione(void)
{
  for(int ix=0;ix<Q;ix++)
    for(int iy=0;iy<Q;iy++)
      for(int i=0; i<Q; i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i] = fnew[ix][iy][i];
}
void LatticeBoltzmann::Inicie(double rho0, double Jx0, double Jy0)
{
  for(int ix=0;ix<Q;ix++)
    for(int iy=0;iy<Q;iy++)
      for(int i=0; i<Q; i++)
	f[ix][iy][i]=feq(rho0, Jx0, Jy0,i);
}
void LatticeBoltzmann::ImponerCampos(int t)
{
  int ix=Lx/2, iy=Ly/2;
  double lambda,omega,rho0,Jx0,Jy0;
  
  lambda=10;
  omega=(2*M_PI*C)/lambda;
  rho0=10*sin(omega*t);
  Jx0=Jx(ix,iy,false);
  Jy0=Jy(ix,iy,false);
  
  for(int i=0;i<Q;i++)
    fnew[Lx/2][Ly/2][i]=feq(rho0,Jx0,Jy0,i);
}
void LatticeBoltzmann::Imprimase(const char * NombreArchivo)
{
  ofstream MiArchivo(NombreArchivo);
  double rho0;       
  for(int ix=0;ix<Lx;ix++)
    {
    for(int iy=0;iy<Ly;iy++)
      {
	rho0=rho(ix,iy,true);
	MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
      }
    MiArchivo<<endl;
    }
  MiArchivo.close();
}
int main(void)
{
  LatticeBoltzmann Ondas;
  int tmax = 100;

  Ondas.Inicie(0,0,0);
  
  for(int t=0;t<tmax;t++)
    {
      Ondas.Colisione();
      Ondas.ImponerCampos(t);
      Ondas.Adveccione();
    }
  Ondas.Imprimase("Ondas.dat");
  
  return 0;
}
