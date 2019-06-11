//Mi Primer Programa en C++
#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;
//---------- Tamaño -----------
const double Lx=60;
const double Ly=120;
const int Nx=5,Ny=5, N=Nx*Ny;
//---------- Constantes de fuerzas --------
const double G=1.0;
const double K=1.0e4;
const double epsilon =1.0;
const double r0 = 5;
//---------- PEFRL --------------------
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
  double GetVx(void){return V.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double GetVy(void){return V.y();}; //Inline
  double Getz(void){return r.z();}; //Inline
  double GetVz(void){return V.z();}; //Inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0){
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
  void CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2);
  void Paredes(Cuerpo & Molecula);
  void CalculeTodasLasFuerzas(Cuerpo * Molecula);
};
void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2){
  vector3D dr=Molecula2.r-Molecula1.r;
  double Norma_dr = norma(dr);
  double r0sobrer = r0/Norma_dr;
  double r0sobrer6 = pow(r0sobrer,6);
  double aux=-12*epsilon*r0sobrer6*(r0sobrer6-1)/(Norma_dr*Norma_dr);
  vector3D F1=dr*aux;
  Molecula1.AgregueFuerza(F1);  Molecula2.AgregueFuerza(F1*(-1));
}
void Colisionador::Paredes(Cuerpo & Molecula){
  double h,R = Molecula.R;
  vector3D x,y,F;
  x.cargue(1,0,0);
  y.cargue(0,1,0);
  double getx = Molecula.Getx();
  double gety = Molecula.Gety();
  if(getx-R<0){
    h = R-getx;
    F = x*K*pow(h,1.5);
    Molecula.AgregueFuerza(F);
  }
  if(getx+R>Lx){
    h = R + getx - Lx;
    F = x*(-K*pow(h,1.5));
    Molecula.AgregueFuerza(F);
  }
  if(gety-R<0){
    h = R-gety;
    F = y*K*pow(h,1.5);
    Molecula.AgregueFuerza(F);
  }
  if(gety+R>Ly){
    h = R + gety - Ly;
    F = y*(-K*pow(h,1.5));
    Molecula.AgregueFuerza(F);
  }
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Molecula){
  int i,j;
  for(i=0;i<(N+4);i++) Molecula[i].BorreFuerza();
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Molecula[i],Molecula[j]);
  for(i=0;i<N;i++) Paredes(Molecula[i]);
}
//------------------ Funciones Globales -----------------
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'pelicula.gif'"<<endl;
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
  cout<<" , "<<Lx/7<<"*t,0";
  cout<<" , "<<Lx/7<<"*t,"<<Ly;
  cout<<" , 0,"<<Ly/7<<"*t";
  cout<<" , "<<Lx<<","<<Ly/7<<"*t";
}
void TermineCuadro(void){
  cout<<endl;
}

int main(void){
  Cuerpo Molecula[N+4];
  Colisionador Newton;
  double t,tdibujo,dt=1.0e-3;
  Crandom ran64(1);
  double KT=10;
  double m0=1,R0=1;
  double y,yprom;
  int i,j;
  
  // InicieAnimacion();
  
  double dx=Lx/(Nx+1);
  double dy=10;//Ly/(Ny+1);
  double Rpared=500*(Lx+Ly);
  double Mpared=500*m0;
  //---------------(x0,y0,z0,Vx0,Vy0,Vz0, m0, R0)
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      
      double x0=(i+1)*dx, y0=(j+1)*dy;
      double V0=sqrt(2*KT/m0), theta=2*M_PI*ran64.r();
      double V0x=V0*cos(theta), V0y=V0*sin(theta);
      
      Molecula[i*Ny+j].Inicie(x0,y0,0,V0x,V0y,0,m0,R0);
    }
  }
  
  // -------- pared derecha----------
  // Molecula[N].Inicie(Lx+Rpared,Ly*0.5,0,0,0,0,Mpared,Rpared);
  // -------- pared izquierda----------
  // Molecula[N+1].Inicie(-Rpared,Ly*0.5,0,0,0,0,Mpared,Rpared);
  // -------- pared superior----------
  // Molecula[N+2].Inicie(Lx*0.5,Ly+Rpared,0,0,0,0,Mpared,Rpared);
  // -------- pared inferior----------
  // Molecula[N+3].Inicie(Lx*0.5,-Rpared,0,0,0,0,Mpared,Rpared);
  
  double T=150;
  for(t=tdibujo=0;t<T;t+=dt,tdibujo+=dt){
	  
	//------------ promedia alturas
    // for(i=0;i<N;i++) y += Molecula[i].Gety();
	// yprom = y / N;
	// y = 0;
	// cout<<t<<" "<<yprom<<endl;
    
    // if(tdibujo>1/24){
      // InicieCuadro();
      // for(i=0;i<N;i++) Molecula[i].Dibujese();
      // TermineCuadro();
      // tdibujo=0;
    // }
	
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
  
  T=200;
  for(t=tdibujo=0;t<T;t+=dt,tdibujo+=dt){
	  
	for(i=0;i<N;i++) cout<<Molecula[i].GetVy()<<endl;
	
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
  
  
  return 0;
}

