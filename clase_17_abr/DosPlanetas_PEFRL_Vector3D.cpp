//Mi Primer Programa en C++
#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
using namespace std;

const int N=2;

const double G=1.0;

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
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0,double R0);
  void BorreFuerza(void){F.cargue(0,0,0);};
  void AgregueFuerza(vector3D dF){F+=dF;};
  void Mueva_r(double dt,double Coeficiente);
  void Mueva_V(double dt,double Coeficiente);
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
void Cuerpo::Mueva_r(double dt,double Coeficiente){
  r+=V*(dt*Coeficiente);
}
void Cuerpo::Mueva_V(double dt,double Coeficiente){
  V+=F*(dt*Coeficiente/m);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//------------------ Clase Colisionador -----------------
class Colisionador{
private:
public:
  void CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2);
  void CalculeTodasLasFuerzas(Cuerpo * Planeta);
};
void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2){
  vector3D dr=Planeta2.r-Planeta1.r;
  double aux=G*Planeta1.m*Planeta2.m*pow(norma2(dr),-1.5);
  vector3D F1=dr*aux;
  Planeta1.AgregueFuerza(F1);  Planeta2.AgregueFuerza(F1*(-1));
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Planeta){
  int i,j;
  for(i=0;i<N;i++) Planeta[i].BorreFuerza();
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Planeta[i],Planeta[j]);
}

//------------------ Funciones Globales -----------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'pelicula.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-150:150]"<<endl;
  cout<<"set yrange[-150:150]"<<endl;
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

int main(void){
  Cuerpo Planeta[N];
  Colisionador Newton;
  double t,tdibujo,dt=1.0;
  double r=100,m0=10,m1=1;
  double M,omega,T,x0,x1,V0,V1;
  int i;

  InicieAnimacion();

  M=m0+m1; omega=sqrt(G*M/(r*r*r)); T=2*M_PI/omega;
  x0=m1*r/M; x1=x0-r; V0=omega*x0; V1=omega*x1;
  //---------------(x0,y0,z0,Vx0,Vy0,Vz0, m0, R0)
  Planeta[0].Inicie(x0, 0, 0,  0, V0,  0, m0, 10);
  Planeta[1].Inicie(x1, 0, 0,  0, V1,  0, m1,  5);

  for(t=tdibujo=0;t<2*T;t+=dt,tdibujo+=dt){
    //cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
    
    if(tdibujo>T/1000){
      InicieCuadro();
      for(i=0;i<N;i++) Planeta[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    //Moverlo Segun PEFRL Orden 4
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Zi);
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Coeficiente1);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Coeficiente2);
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Coeficiente1);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Zi);
  }
  
  return 0;
}
