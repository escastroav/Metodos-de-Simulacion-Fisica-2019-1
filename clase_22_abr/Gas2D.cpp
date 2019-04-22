//Mi Primer Programa en C++
#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;

const double K = 1e4;

//Dimensiones de las paredes
const double Lx = 60, Ly = 60;

//Parametro de arreglo de granos en pos. inicial
const int Nx = 5, Ny = 5;
const int N = Nx*Ny;

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
  void CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2);
  void CalculeTodasLasFuerzas(Cuerpo * Molecula);
};
void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2){
  vector3D dr=Molecula2.r-Molecula1.r;
  double Norma_dr = norma(dr);
  //vector3D F2;F2.cargue(0,0,0);
  double s = (Molecula1.R + Molecula2.R) - Norma_dr;
  if(s > 0) //hay choque
    {
      vector3D F2 = dr*(K*pow(s,1.5)/Norma_dr); 
      Molecula2.AgregueFuerza(F2);  Molecula1.AgregueFuerza(F2*(-1));
    }
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Molecula){
  int i,j;
  for(i=0;i<(N+4);i++) Molecula[i].BorreFuerza(); //Incluyen paredes en N+4
  for(i=0;i<N;i++)
    for(j=i+1;j<(N+4);j++) //Incluyen paredes QUE RECIBEN la fuerza
      CalculeFuerzaEntre(Molecula[i],Molecula[j]);
}

//------------------ Funciones Globales -----------------
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl; 
  //cout<<"set output 'pelicula.gif'"<<endl;
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
    //Dibujando paredes
    cout << " , " << Lx/7 << "*t,0" //pared abajo
	 << " , " << Lx/7 << "*t,"<<Ly //pared arriba
	 <<" , 0,"<< Ly/7 << "*t" //pared izquierda
	 <<" , " << Lx << " , " <<Ly/7<<"*t"; //pared derecha
}
void TermineCuadro(void){
    cout<<endl;
}

int main(void){
  Cuerpo Molecula[N+4];
  Colisionador Newton;
  Crandom ran64(1);
  double t,tdibujo,dt=0.1;
  int i,j;
  double m0 = 1;
  double R0 = 2;
  double kT = 0.01; //K_b * T definido por energia cinetica promedio


  InicieAnimacion();
  
  //Definir Paredes
  double Mpared = 100*m0, Rpared = 1e4;
  //----------------(x0,y0         ,z0,Vx0,Vy0,Vz0, m0, R0)
  //Pared arriba
  Molecula[N].Inicie(Lx/2,Ly+Rpared, 0,0  ,0  ,0  ,Mpared,Rpared);
  //Pared abajo
  Molecula[N+1].Inicie(Lx/2,-Rpared, 0,0  ,0  ,0  ,Mpared,Rpared);
  //Pared derecha
  Molecula[N+2].Inicie(Lx+Rpared,Ly/2, 0,0  ,0  ,0  ,Mpared,Rpared);
  //Pared izquierda
  Molecula[N+3].Inicie(-Rpared,Ly/2, 0,0  ,0  ,0  ,Mpared,Rpared);


  //Definir pos inicial de granos
  double dx = Lx/(Nx+1), dy = Ly/(Ny+1), x0, y0, Vx0, Vy0;
  double V0 = sqrt(2*kT/m0); double theta;
  for(i=0;i<Nx;i++)
    for(j=0;j<Ny;j++)
      {
	x0 = (i+1)*dx;
	y0 = (j+1)*dy;
	theta = 2*M_PI*ran64.r();
	Vx0 = V0*cos(theta);
	Vy0 = V0*sin(theta);
	//---------------------(x0,y0,z0,Vx0,Vy0,Vz0, m0, R0)
	Molecula[i*Ny+j].Inicie(x0,y0,0,Vx0,Vy0, 0, m0, R0);
      }
  
  

  
  //Correr simulacion
  double T = 200;
  for(t=tdibujo=0;t<T;t+=dt,tdibujo+=dt){
    //for(i=0;i<N;i++)cout<<Molecula[i].Getx()<<" "<<Molecula[i].Gety()<<endl;
    
    if(tdibujo>T/1000){
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
  
  //Imprimir resultados
      
  
  return 0;
}
