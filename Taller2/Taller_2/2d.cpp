#include<iostream>
#include<cstdlib>
#include<cmath>
#include<fstream>

const int D=2;
const int Q = 5;

const int Lx=400;
const int Ly=200;

const double W0=1.0/3;//Weight for the zero vector

const double C=0.5;//Velocidad de las ondas c<1/sqrt(2) para que la funci�n de equilibrio no de cero.

const double tau=0.5;//Relaxation time
const double Utau = 1.0/tau;
const double UmUtau=1-Utau;

class LatticeBoltzmann{
private:
  double w[Q];
  int v[D][Q];//V[0][i]=V_ix, V[1][i]=V_iy
  double f[Lx][Ly][Q]; double fnew[Lx][Ly][Q];
public:
  LatticeBoltzmann(void);
  double CCelda(int ix, int iy);
  double rho(int ix, int iy, bool UsarNew);
  double Jx(int ix, int iy, bool UsarNew);
  double Jy(int ix, int iy, bool UsarNew);
  double feq(int i, double rho0, double Jx0, double Jy0, double Ci);//we don't use double rho so we don't confuse the value with the function rho()
  void Inicie(double rho0, double Jx0, double Jy0);
  void Colisione(void);
  void Adveccione(void);
  void ImponerCampos(int t);
  void Imprimase(const char* NombreDelArchivo);
};

LatticeBoltzmann::LatticeBoltzmann(void){//Funci�n Constructor de la clase
  //Pesos
  w[0]=W0;
  w[1]=w[2]=w[3]=w[4]=(1-W0)/4.0;
  //Velocidades
  v[0][0]=0;//v_0x
  v[1][0]=0;//v_0y
  v[0][1]=1; v[1][1]=0;
  v[0][2]=0; v[1][2]=1;
  v[0][3]=-1; v[1][3]=0;
  v[0][4]=0; v[1][4]=-1;
}

double LatticeBoltzmann::CCelda(int ix, int iy){
  //  return 0.5;
  return 0.25*C*(tanh(100-ix)+3);
}

double LatticeBoltzmann::rho(int ix, int iy, bool UsarNew){
  int i;
  double suma=0.0;
  if(UsarNew)
    for(i=0;i<Q;i++)
      suma+=fnew[ix][iy][i];
  else
    for(i=0;i<Q;i++)
      suma+=f[ix][iy][i];
  return suma;
}

double LatticeBoltzmann::Jx(int ix, int iy, bool UsarNew){
  int i;
  double suma=0.0;
  if(UsarNew)
    for(i=0;i<Q;i++)
      suma+=v[0][i]*fnew[ix][iy][i];
  else
    for(i=0;i<Q;i++)
      suma+=v[0][i]*f[ix][iy][i];
  return suma;
}

double LatticeBoltzmann::Jy(int ix, int iy, bool UsarNew){
  int i;
  double suma=0.0;
  if(UsarNew)
    for(i=0;i<Q;i++)
      suma+=v[1][i]*fnew[ix][iy][i];
  else
    for(i=0;i<Q;i++)
      suma+=v[1][i]*f[ix][iy][i];
  return suma;
}

double LatticeBoltzmann::feq(int i, double rho0, double Jx0, double Jy0, double Ci){
  double TresCi2 = 3*Ci*Ci;
  double aux=1-TresCi2*(1-W0);
  if(i==0)
    return rho0*aux;
  else
    return w[i]*(TresCi2*rho0+3*(v[0][i]*Jx0+v[1][i]*Jy0));
}

void LatticeBoltzmann::Inicie(double rho0, double Jx0, double Jy0){
  double Ci;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      Ci=CCelda(ix,iy);
      for(int i=0;i<Q;i++)
	f[ix][iy][i]=fnew[ix][iy][i]=feq(i,rho0,Jx0,Jy0,Ci);
    }
}

void LatticeBoltzmann::Colisione(void){
  double rho0, Jx0, Jy0, Ci;
  //Para cada celda
  for(int ix=0; ix<Lx;ix++)
    for(int iy=0; iy<Ly;iy++){
      //----Calculamos las cantidades macrosc�picas
      rho0=rho(ix,iy, false); Jx0=Jx(ix,iy, false); Jy0=Jy(ix,iy, false); Ci=CCelda(ix,iy);
      for(int i=0;i<Q;i++)
	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(i,rho0,Jx0,Jy0,Ci);
    }
}

void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[(ix+v[0][i]+Lx)%Lx][(iy+v[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];//Condiciones de frontera peri�dicas
}

void LatticeBoltzmann::ImponerCampos(int t){//Colocamos un palito que sube y baja en el centro del lattice
  double A = 10.0, omega, lambda=10.0;//amplitude, frecuency and wavelength
  omega = 2*M_PI*C/lambda;
  double rho0, Jx0, Jy0;
  int i, ix=0, iy;
  //A plane wave goes like A*exp(i(vec(k)*vec(x)-omega*t)), since x=0, we use sin
  rho0=A*sin(omega*t);//We need a function that starts at 0 so it's a smooth transition
  for(iy=0;iy<Ly;iy++){
    Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
    for(i=0;i<Q;i++)
      fnew[ix][iy][i]=feq(i, rho0, Jx0, Jy0, C);
  }
}

void LatticeBoltzmann::Imprimase(const char* NombreDelArchivo){
  std::ofstream File(NombreDelArchivo);
  double rho0;
  for(int ix=0; ix<201; ix++){
    for(int iy=0; iy<201; iy++){
      rho0=rho(ix,iy,true);
      File<<ix<<" "<<iy<<" "<<rho0<<"\n";
    }
    File<<std::endl;
  }
  File.close();
}

int main(){
  int t, tmax=400;
  double RHO0=0.0;double JX0=0;double JY0=0;
  
  LatticeBoltzmann Ondas;
  Ondas.Inicie(RHO0,JX0,JY0);

  for(t=0;t<tmax;t++){
    Ondas.Colisione();
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
  }

  Ondas.Imprimase("2d.dat");

  return 0;
}
