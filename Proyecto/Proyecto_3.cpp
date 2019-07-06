#include<iostream>
#include<cmath>
#include<vector>
#include "Vector.h"
#include "Random64.h"


const double g = 9.8, k = 1e4, Gamma = 5e2, kcundall = 1e308, MU = 0.4;
const double lx = 0, ly = 110;

const int N = 402, NumP = 430;

const double E = 0.1786178958448091;
const double L = -0.2123418310626054;
const double X = -0.06626458266981849;
const double L1 = (1-2*L)*0.5;
const double XE = 1-2*(X*E);
//----------- Delcaracion de las clases-------

class Cuerpo;
class Colisionador;

//----------- Clase cuerpo -------------------

class Cuerpo{
  //Los datos en cuerpo son privados, las ordenes son publicas
private:     
  double m,R,theta,I;
  vector3D r, V, F, omega, tau;

public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0, double Vz0,
	      double theta0, double omega0,
	      double m0,double R0);
  void BorreFuerzaYTorque(void){F.cargue(0,0,0);tau.cargue(0,0,0);}
  void AgregueFuerza(vector3D dF){F+=dF;}
  void AgregueTorque(vector3D dtau){tau+=dtau;}
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_v(double dt, double Coeficiente);
  void Dibujese(void);

  //Se hacen funciones separadas macro, i.e, que vaya haciendo la funcion mientras se compila el programa, para sacar los datos privados y solo imprimirlos y no modificarlos
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();}; //Inline
  double GetVx(void){return V.x();};
  friend class Colisionador;
};
//---------------- Funciones de la calse Cuerpo----------------

//Para que el compilador sepa que una funcion es de la clase cuerpo, se le escribe lo del principio, y puede usar las variables de la clase cuerpo

//Se definen las condiciones iniciales del cuerpo
void Cuerpo::Inicie(double x0,double y0,double z0,
		    double Vx0,double Vy0, double Vz0,
		    double theta0, double omega0,
		    double m0,double R0){
  m=m0; R=R0; theta=theta0; I=(2.0/5)*m*R*R;
  r.cargue(x0,y0,z0);
  V.cargue(Vx0,Vy0,Vz0);
  omega.cargue(0,0,omega0);
}

//Se mueve r y v
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(dt*Coeficiente); theta+=omega.z()*dt*Coeficiente; 
}
void Cuerpo::Mueva_v(double dt, double Coeficiente){
  V+=F*(dt*Coeficiente/m); omega+=tau*(dt*Coeficiente/I);
}

void Cuerpo::Dibujese(void){
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

//-------------- Clase Colisionador ---------

class Colisionador{
private:
  vector3D ele[N][N]; bool EstabaEnContacto[N][N];
  vector3D elePared[N][NumP]; bool EstabaEnContactoPared[N][NumP];
public:
  void Inicie(void);
  void CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2, vector3D & ele, bool & EstabaEnContacto, double dt);
  void CalculeTodasLasFuerzas(std::vector<Cuerpo> & Molecula, std::vector<Cuerpo> & Paredes, double dt);
};

//---------------- Funciones de la calse Cuerpo----------------

void Colisionador::Inicie(void){
  for(int i = 0 ; i < N ; i++){
    for(int j = 0; j < N; j++){
      ele[i][j].cargue(0,0,0);
      EstabaEnContacto[i][j]=false;
    }
    for(int k = 0; k < NumP; k++){
      elePared[i][k].cargue(0,0,0);
      EstabaEnContactoPared[i][k]=false;
    }
  }
}

void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2, vector3D & ele, bool & EstabaEnContacto, double dt){
  double ERRF = 1e-8;
  vector3D dr = Molecula2.r - Molecula1.r; double norma_dr = norma(dr);
  double s = (Molecula2.R + Molecula1.R) - norma_dr;
  if(s>0){
    vector3D n = dr/norma_dr, Vc = Molecula2.V-Molecula1.V-((Molecula1.omega*Molecula1.R + Molecula2.omega*Molecula2.R)^n);
    double Componente_Vcn = Vc*n;
    vector3D Vcn = n*Componente_Vcn, Vct = Vc - Vcn;
    double Norma_Vct = norma(Vct);
    vector3D t;
    if(Norma_Vct<ERRF) t.cargue(0,0,0); else t = Vct/Norma_Vct;
    double ComponenteFn;
    ComponenteFn = k*pow(s,1.5);
    double m12 = Molecula1.m*Molecula2.m/(Molecula1.m+Molecula2.m);
    ComponenteFn-=m12*sqrt(s)*Gamma*Componente_Vcn; if(ComponenteFn < 0) ComponenteFn = 0;
    vector3D Fn = n*ComponenteFn;

    ele+=Vct*dt; vector3D Ft = ele*(-kcundall);
    double FtMax = MU*std::fabs(ComponenteFn);
    if(norma(Ft)>FtMax){
      Ft = ele*(-FtMax/norma(ele));
    }
    
    vector3D F2 = Fn + Ft;
    vector3D T2 = ((-Molecula2.R)*n)^F2;
    vector3D T1 = ((-Molecula1.R)*n)^F2;
    Molecula2.AgregueFuerza(F2); Molecula1.AgregueFuerza(F2*(-1));
    Molecula2.AgregueTorque(T2); Molecula1.AgregueTorque(T1);

    EstabaEnContacto = true;
  }

  else if(EstabaEnContacto){
    ele.cargue(0,0,0);
    EstabaEnContacto=false;
  }
}

void Colisionador::CalculeTodasLasFuerzas(std::vector<Cuerpo> & Molecula, std::vector<Cuerpo> & Paredes, double dt){
  int i,j,k; vector3D Gvector;
  
  for(i=0;i<Molecula.size();i++){
    Molecula[i].BorreFuerzaYTorque();
  }

  for(i=0;i<Paredes.size();i++){
    Paredes[i].BorreFuerzaYTorque();
  }

  Gvector.cargue(0,-g,0);
  for(i=0;i<Molecula.size();i++){
    Molecula[i].AgregueFuerza(Gvector*Molecula[i].m);
  }
  
  for(i=0;i<Molecula.size();i++){
    for(k=0;k<Paredes.size();k++){
      CalculeFuerzaEntre(Molecula[i],Paredes[k],elePared[i][k],EstabaEnContactoPared[i][k],dt);
    }
    for(j=i+1;j<Molecula.size();j++){
      CalculeFuerzaEntre(Molecula[i],Molecula[j],ele[i][j],EstabaEnContacto[i][j],dt);
    }
  }
}

//---------------- Funciones Globales ----------------

void InicieAnimacion(void){
  //std::cout<<"set terminal gif animate"<<std::endl;
  //std::cout<<"set output 'pelicula.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange[-100:"<<lx+100<<"]"<<std::endl;
  std::cout<<"set yrange[-10:"<<ly<<"]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange[0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}

void InicieCuadro(void){
  std::cout<<"plot 0,0 ";
  std::cout<<" , "<<(lx+172)/7<<"*t-86.0,0";
  std::cout<<" , "<<(lx+172)/7<<"*t-86.0,101";
}

void TermineCuadro(void){
  std::cout<<std::endl;
}


int main(){
  
  //Se define una variable de cuerpo llamada Molecula
  std::vector<Cuerpo> Molecula(1);
  std::vector<Cuerpo> Pared(NumP);
  Colisionador Newton;
  Crandom ran64(1);
  double t,tdibujo;
  int tcrear;
  double dt=0.001;
  int i,j;
  double m0 = 1.0, R0 = 1.0;
  double kT = 10.0;
  double V0 = sqrt(2*kT/m0), deltax;
  double h = 100;
    
  InicieAnimacion();

  //---------- Primera molecula -----------
  deltax = 0.02*ran64.r() - 0.01;
  Molecula[0].Inicie(lx/2.0 +deltax,h,0,0,0,0,0,0,m0,R0); 

  //---------- Paredes ---------------
  double Mpared=100*m0, Rpared = 0.2;
  Cuerpo MolPared;
  for(int i = 0; i < NumP; i++){
    Pared[i].Inicie(-86 + 2*i*Rpared,-Rpared,0,0,0,0,0,0,Mpared,Rpared);
  }

  Cuerpo NuevaMol; 
  double T = 401;
    //Se empieza a mover el cuerpo
  for(t=tdibujo=0,tcrear=1;t<T;t+=dt,tdibujo+=dt,tcrear++){

    //----crear una nueva movecula cada cierto tiempo-------
    if(tcrear == 1000){
      deltax = 0.02*ran64.r() - 0.01;
      NuevaMol.Inicie(lx/2.0 +deltax,h,0,0,0,0,0,0,m0,R0);  
      Molecula.push_back(NuevaMol);
      tcrear = 1;
    }
    
    if(tdibujo>T/800){ //800 frames
      InicieCuadro();
      for(i=0;i<Molecula.size();i++){
	Molecula[i].Dibujese();
      }
      TermineCuadro();
      tdibujo=0;
    }
    
    //Mover segunn F-R orden 4
    for(i=0;i<Molecula.size();i++) Molecula[i].Mueva_r(dt,E);
    
    Newton.CalculeTodasLasFuerzas(Molecula,Pared,dt);
    for(i=0;i<Molecula.size();i++) Molecula[i].Mueva_v(dt,L1);
    
    for(i=0;i<Molecula.size();i++) Molecula[i].Mueva_r(dt,X);
    
    Newton.CalculeTodasLasFuerzas(Molecula,Pared,dt);
    for(i=0;i<Molecula.size();i++) Molecula[i].Mueva_v(dt,L);
    
    for(i=0;i<Molecula.size();i++) Molecula[i].Mueva_r(dt,XE);
    
    Newton.CalculeTodasLasFuerzas(Molecula,Pared,dt);
    for(i=0;i<Molecula.size();i++) Molecula[i].Mueva_v(dt,L);
    
    for(i=0;i<Molecula.size();i++) Molecula[i].Mueva_r(dt,X);

    Newton.CalculeTodasLasFuerzas(Molecula,Pared,dt);
    for(i=0;i<Molecula.size();i++) Molecula[i].Mueva_v(dt,L1);

    for(i=0;i<Molecula.size();i++) Molecula[i].Mueva_r(dt,E);
  }

  //------ imprime resultado -----------------------

  /*for(int n=0;n<N;n++)
    std::cout<<Molecula[n].GetVx()<<std::endl;*/
  return 0;
}
