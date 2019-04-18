#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"

using namespace std;

const int N = 2;

const double G=1.0;

const double Zi = 0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Xi = -0.06626458266981849;

const double Coeff1 = (1-2*Lambda)*0.5;
const double Coeff2 = (1-2*(Xi+Zi));

//------Declaracion de clases
class Cuerpo;
class Colisionador;

class Cuerpo{
private:
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0,double R0);
  void BorreFuerza(void){F.cargue(0,0,0);};
  void AgregueFuerza(vector3D dF){F+=dF;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_v(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();}; //Inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0,double R0){
  r.cargue(x0,y0,z0); V.cargue(Vx0,Vy0,Vz0); m=m0; R=R0;
}

void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(dt*Coeficiente);
}
void Cuerpo::Mueva_v(double dt, double Coeficiente){
  V+=F*(dt*Coeficiente/m);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//------------------ Clase Colisionador -----------------
class Colisionador
{
private:
  
public:
  void CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2);
  void CalculeTodasLasFuerzas(Cuerpo * Planeta);
};
void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2)
{
  vector3D dr=Planeta2.r-Planeta1.r;
  double aux = G*Planeta1.m*Planeta2.m*pow(norma2(dr),-1.5);
  vector3D F1 = dr*aux;
  Planeta1.AgregueFuerza(F1);
  Planeta2.AgregueFuerza(F1*(-1));
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Planeta)
{
  for(int i=0; i < N; i++)
    Planeta[i].BorreFuerza();
  for(int i=0; i < N; i++)
    for(int j=i+1; j < N; j++)
      CalculeFuerzaEntre(Planeta[i],Planeta[j]);
}
//------------------ Funciones Globales -----------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'pelicula.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-1500:1500]"<<endl;
  cout<<"set yrange[-1500:1500]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}

void PEFRL (double dt, Cuerpo * Planeta, Colisionador & Newton)
{
  for(int i = 0; i < N; i++) Planeta[i].Mueva_r(dt,Zi);
  Newton.CalculeTodasLasFuerzas(Planeta);
  for(int i = 0; i < N; i++) Planeta[i].Mueva_v(dt,Coeff1);
  for(int i = 0; i < N; i++) Planeta[i].Mueva_r(dt,Xi);
  Newton.CalculeTodasLasFuerzas(Planeta);
  for(int i = 0; i < N; i++) Planeta[i].Mueva_v(dt,Lambda);
  for(int i = 0; i < N; i++) Planeta[i].Mueva_r(dt,Coeff2);    
  Newton.CalculeTodasLasFuerzas(Planeta);    
  for(int i = 0; i < N; i++) Planeta[i].Mueva_v(dt,Lambda);    
  for(int i = 0; i < N; i++) Planeta[i].Mueva_r(dt,Xi);    
  Newton.CalculeTodasLasFuerzas(Planeta);    
  for(int i = 0; i < N; i++)Planeta[i].Mueva_v(dt,Coeff1);    
  for(int i = 0; i < N; i++) Planeta[i].Mueva_r(dt,Zi);
}

int main(void){
  Cuerpo Planeta[N];
  Colisionador Newton;
  double t,tdibujo,dt=0.1;
  double r=1000,
    m0 = 1047,
    m1 = 1,
    M,
    x0,
    x1,
    omega,
    T,
    V0,
    V1;

  //InicieAnimacion();

  M = m0+m1;
  omega = sqrt(G*M/(r*r*r));
  T = 2*M_PI/omega;
  
  x0 = m1*r/M;
  x1 = x0-r;
  V0 = omega*x0;
  V1 = omega*x1;
  
  //------------(x0,y0,z0,Vx0,Vy0,Vz0,  m0, R0)
  Planeta[0].Inicie(x0, 0, 0, 0, V0, 0, m0, 100);
  Planeta[1].Inicie(x1, 0, 0, 0, V1, 0, m1, 50);
  
  for(t=tdibujo=0;t<20*T;t+=dt,tdibujo+=dt){
    for(int i = 0; i < N; i++)
      cout<<Planeta[i].Getx()
	  <<"\t"
	  <<Planeta[i].Gety()<<"\t";
    cout<<" "<<endl;
    /*if(tdibujo>T/1000){
      InicieCuadro();
      for(int i = 0; i < N; i++) Planeta[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
      }*/
    //Moverse segun PEFRL
    PEFRL(dt,Planeta,Newton);
  }
  
  return 0;
}
