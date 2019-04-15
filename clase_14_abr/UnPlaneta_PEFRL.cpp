#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"

using namespace std;

const double GM=1.0;

const double Theta = 1.0/(2.0-pow(2.0,1.0/3.0));
const double ThetaMed = Theta/2;
const double UnoThetaMed = 1-ThetaMed;
const double UnoDosTheta = 1 - 2*Theta;

class Cuerpo{
private:
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0,double R0);
  void CalculeFuerza(void);
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_v(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();}; //Inline
};
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0,double R0){
  r.cargue(x0,y0,z0); V.cargue(Vx0,Vy0,Vz0); m=m0; R=R0;
}
void Cuerpo::CalculeFuerza(void){
  double aux=-GM*m*pow(norma2(r),-1.5);
  F=r*aux;

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
  Cuerpo Planeta;
  double t,tdibujo,dt=1.0;
  double r=100,omega,T,V;

  InicieAnimacion();

  omega=sqrt(GM/(r*r*r)); T=2*M_PI/omega; V=omega*r;
  //------------(x0,y0,z0,Vx0,Vy0  ,Vz0,  m0, R0)
  Planeta.Inicie(r ,0 , 0,  0.1*V,0.5*V,  0, 1.0, 5);
  
  for(t=tdibujo=0;t<1.1*T;t+=dt,tdibujo+=dt){
    //cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
    
    if(tdibujo>T/10000){
      InicieCuadro();
      Planeta.Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    //Moverse segun ForestRuth
    Planeta.Mueva_r(dt,ThetaMed);
    
    Planeta.CalculeFuerza();
    Planeta.Mueva_v(dt,Theta);
    
    Planeta.Mueva_r(dt,UnoThetaMed);
    
    Planeta.CalculeFuerza();
    Planeta.Mueva_v(dt,UnoDosTheta);
    
    Planeta.Mueva_r(dt,UnoThetaMed);
    
    Planeta.CalculeFuerza();
    Planeta.Mueva_v(dt,Theta);

    Planeta.Mueva_r(dt,ThetaMed);
  }
  
  return 0;
}
