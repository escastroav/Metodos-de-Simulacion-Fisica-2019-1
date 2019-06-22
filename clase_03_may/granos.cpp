//Mi Primer Programa en C++
#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;

const double g=9.81, K=1.0e4, Gamma=50, Kcundall=10,MU=0.4;
const double Lx=60,Ly=60;
const int Nx=5,Ny=5, N=Nx*Ny;

const double Zi=0.1786178958448091e0;
const double Lambda=0.2123418310626054*(-1);
const double Xi=0.06626458266981849*(-1);

const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=(1-2*(Xi+Zi));

//------------ Declaracion de las clases --------
class Cuerpo;
class Colisionador;
//------------ Clase Cuerpo --------
class Cuerpo{
private:
  vector3D r,V,F,omega,tau; double m,R, theta, I;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double theta0, double omega0, double m0, double R0);
  void BorreFuerzaYTorque(void){F.cargue(0,0,0);tau.cargue(0,0,0);};
  void AgregueFuerza(vector3D dF){F+=dF;};
  void AgregueTorque(vector3D dtau){tau+=dtau;};
  void Mueva_r(double dt,double Coeficiente);
  void Mueva_V(double dt,double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();}; //Inline
  double GetVx(void){return V.x();}; //Inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double theta0, double omega0, double m0, double R0){
  r.cargue(x0,y0,z0); V.cargue(Vx0,Vy0,Vz0); omega.cargue(0,0,omega0);
  theta = theta0; m=m0; R=R0; I=2.0/5*0*R*R;
}
void Cuerpo::Mueva_r(double dt,double Coeficiente){
  r+=V*(dt*Coeficiente);
  theta+=omega.z()*(dt*Coeficiente);
}
void Cuerpo::Mueva_V(double dt,double Coeficiente){
  V+=F*(dt*Coeficiente/m);
  omega+=tau*(dt*Coeficiente/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}
//------------------ Clase Colisionador -----------------
class Colisionador{
private:
  vector3D ele[N+4][N+4]; bool EstabaEnContacto[N+4][N+4];
public:
  void Inicie(void);
  void CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2);
  void CalculeTodasLasFuerzas(Cuerpo * Molecula);
};
void Colisionador::Inicie(void)
{
  for(int i=0;i<N+4;i++)
    for(int j=0;j<N+4;j++)
      {
	ele[i][j].cargue(0,0,0);
	EstabaEnContacto[i][j]=false;
      }
    
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2, vector3D & ele, bool EstabaEnContacto){
  double ERRF = 1e-8;
  vector3D dr=Molecula2.r-Molecula1.r; double Norma_dr=norma(dr);
  double s=(Molecula2.R+Molecula1.R)-Norma_dr;
  if(s>0)
    { //Si hay choque
      //Calcular cantidades auxiliares
      vector3D
	n = dr/Norma_dr,
	Vc = Molecula2.V-Molecula1.V
	-((Molecula1.omega*Molecula1.R+Molecula2.omega*Molecula2.R)^n),
	Vcn = (Vc*n)*n,
	Vct = Vc - Vcn;
      double Norma_Vct=norma(Vct);
      vector3D t;
      if(Norma_Vct < ERRF)
	t.cargue(0,0,0);
      else
	t=Vct/Norma_Vct;
      //FUERZAS NORMALES
      double ComponenteFn;
      //Fuerza de Hertz
      ComponenteFn=K*pow(s,1.5);
      //Fuerza disipativa
      double m12=(Molecula1.m*Molecula2.m)/(Molecula1.m+Molecula2.m);
      ComponenteFn-=m12*sqrt(s)*Gamma*Componente_Vcn; if(ComponenteFn<0) ComponenteFn=0;
      //Fuerza normal total
      vector3D Fn=n*ComponenteFn;
      //FUERZAS TANGENCIALES
      //Fuerza de friccion estatica
      ele+=Vct*dt;
      vector3D Ft=ele*(-Kcundall);
      //Fuerza de friccion cinetica
      double FtMax=MU*fabs(ComponenteFn);
      if(norma(Ft)>FtMax) Ft=ele*(-FtMax/norma(ele));
      
      //ENSAMBLAR LA FUERZA TOTAL Y CARGARLA
      vector3D F2 = Fn+Ft;
      Molecula2.AgregueFuerza(F2); Molecula1.AgregueFuerza((-1)*F2);
      vector3D nCruzF2 = n^F2;
      Molecula2.AguregueTorque(nCruzF2*(-Molecula2.R));
      Molecula1.AguregueTorque(nCruzF2*(-Molecula1.R));

      EstabaEnContacto=true;
    }
  else if(EstabaEnContacto)
    {
      ele.cargue(0,0,0);
      EstebaEnContacto=false;
    }
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Molecula){
  int i,j;
  vector3D g_vector;g_vector.cargue(0,-g,0);
  for(i=0;i<(N+4);i++) Molecula[i].BorreFuerzaYTorque();
  for(i=0;i<N;i++)
    for(j=i+1;j<(N+4);j++)
      CalculeFuerzaEntre(Molecula[i],Molecula[j],ele[i][j],EstabaEnContacto[i][j]);
}

//------------------ Funciones Globales -----------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'pelicula.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}

int main(void){
  Cuerpo Molecula[N+4];
  Colisionador Newton;
  Crandom ran64(1);
  double t,tdibujo,dt=1.0e-3;
  double m0=1,R0=2;
  int i,j,n;
  double kT=10;
  

  InicieAnimacion();

  //--------------INICIALIZACION----------------------------

  //PAREDES
  double Mpared=100*m0, Rpared=10000;
  //---------------(x0  ,       y0,z0,Vx0,Vy0,Vz0,m0,R0);
  //Pared arriba
  Molecula[N  ].Inicie(Lx/2,Ly+Rpared, 0, 0, 0, 0, Mpared, Rpared);
  //Pared abajo
  Molecula[N+1].Inicie(Lx/2,  -Rpared, 0, 0, 0, 0, Mpared, Rpared);
  //Pared derecha
  Molecula[N+2].Inicie(Lx+Rpared,Ly/2, 0, 0, 0, 0, Mpared, Rpared);
  //Pared izquierda
  Molecula[N+3].Inicie(  -Rpared,Ly/2, 0, 0, 0, 0, Mpared, Rpared);

  //MOLECULAS
  double dx=Lx/(Nx+1), dy=Ly/(Ny+1), x0,y0;
  double V0=sqrt(2*kT/m0), theta,Vx0,Vy0;
  for(i=0;i<Nx;i++)
    for(j=0;j<Ny;j++){
      x0=(i+1)*dx; y0=(j+1)*dy; 
      theta=2*M_PI*ran64.r(); Vx0=V0*cos(theta); Vy0=V0*sin(theta);
      //---------------------(x0,y0,z0,Vx0,Vy0,Vz0, m0, R0)
      Molecula[i*Ny+j].Inicie(x0,y0, 0,Vx0,Vy0,  0, m0, R0);
    }
      
  //------------CORRER LA SIMULACION-------------------------
  double T=50,Teq=20;
  for(t=tdibujo=0;t<Teq+T;t+=dt,tdibujo+=dt){
    
    if(tdibujo>T/1000){
      // if(t>Teq) for(int n=0;n<N;n++) cout<<Molecula[n].GetVx()<<endl;
      
      InicieCuadro();
      for(i=0;i<N;i++) Molecula[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }
    
    //Moverlo Segun PEFRL Orden 4
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Zi);
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Coeficiente1);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Coeficiente2);
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Coeficiente1);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Zi);
    
  }
  //------------IMPRIMIR RESULTADOS-------------------------
  return 0;
}
