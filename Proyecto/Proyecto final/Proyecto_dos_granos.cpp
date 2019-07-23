#include<iostream>
#include<cmath>
#include<vector>
#include "Vector.h"
#include "Random64.h"


const char FALSE = 'f';
const char TRUE = 't';

const double g = 9.8, k = 1e4, Gamma = 5e2, kcundall = 1e308, MU = 0.4, surfTension=27.5;//surfTension es el gamma de la fuerza capilar
const double lx = 0, ly = 110;

const int N = 50, NumP = 1;
int NN = N*N, NNumP = N*NumP;

const double E = 0.1786178958448091;
const double L = -0.2123418310626054;
const double X = -0.06626458266981849;
const double L1 = (1-2*L)*0.5;
const double XE = 1-2*(X*E);

//const double VaguaMax=0.105;//Volumen de agua máximo entre las partículas Thetamx=30
const double VaguaMax=0.0525;//Volumen de agua máximo entre las partículas Thetamx=25
const double DosVsPi=2*VaguaMax/M_PI;

double Theta2(double s){
  double aux= 1- 2*s/3.0;
  double Theta2=(-s+sqrt(s*s+aux*DosVsPi))/aux;//Usando R=1
  return Theta2;
}

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
  std::vector<vector3D> ele; std::vector<char> EstabaEnContacto;
  std::vector<vector3D> elePared; std::vector<char> EstabaEnContactoPared;
public:
  void Inicie(void);
  void CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2, vector3D & ele, char & Contacto, double dt);
  void CalculeTodasLasFuerzas(std::vector<Cuerpo> & Molecula, std::vector<Cuerpo> & Paredes, double dt);
};

//---------------- Funciones de la calse Cuerpo----------------

void Colisionador::Inicie(void){
  
  ele.resize(NN);
  EstabaEnContacto.resize(NN);
  elePared.resize(NNumP);
  EstabaEnContactoPared.resize(NNumP);
  
  for(int i = 0 ; i < N ; i++){
    for(int j = 0; j < N; j++){
      ele[i*N + j].cargue(0,0,0);
      EstabaEnContacto[i*N + j]=FALSE;
    }
    for(int k = 0; k < NumP; k++){
      elePared[i*NumP + k].cargue(0,0,0);
      EstabaEnContactoPared[i*NumP + k]=FALSE;
    }
  }
}

void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2, vector3D & ele, char & Contacto , double dt){
  double ERRF = 1e-8;
  vector3D dr = Molecula2.r - Molecula1.r; double norma_dr = norma(dr);
  double s = (Molecula2.R + Molecula1.R) - norma_dr;
  if(-s<0.2){//Fuerza capilar
    vector3D Ft, Fn;
    Ft.cargue(0,0,0); Fn.cargue(0,0,0);
    double ComponenteFn = 0;
    vector3D n = dr/norma_dr;
    if(s>0){//Fuerzas de contacto
      vector3D Vc = Molecula2.V-Molecula1.V-((Molecula1.omega*Molecula1.R + Molecula2.omega*Molecula2.R)^n);
      double Componente_Vcn = Vc*n;
      vector3D Vcn = n*Componente_Vcn, Vct = Vc - Vcn;
      double Norma_Vct = norma(Vct);
      vector3D t;
      if(Norma_Vct<ERRF) t.cargue(0,0,0); else t = Vct/Norma_Vct;
      //Hertz
      ComponenteFn += k*pow(s,1.5);
      //Disipación
      double m12 = Molecula1.m*Molecula2.m/(Molecula1.m+Molecula2.m);
      ComponenteFn-=m12*sqrt(s)*Gamma*Componente_Vcn;
      if(ComponenteFn < 0) ComponenteFn = 0;
      //Fricción por Kundall
      ele+=Vct*dt; Ft = ele*(-kcundall);
      double FtMax = MU*std::fabs(ComponenteFn);
      if(norma(Ft)>FtMax){
	Ft = ele*(-FtMax/norma(ele));
      }
      Contacto = TRUE;
    }
    double R = (Molecula1.R+Molecula2.R)/2.0;
    if(s>0) s=0;
    double Th2 =Theta2(-s);
    //Fuerza Capilar
    ComponenteFn-=2*M_PI*R*R*surfTension*Th2/abs(-s-R*Th2);
    Fn = n*ComponenteFn;
    vector3D F2 = Fn + Ft;
    vector3D T2 = ((-Molecula2.R)*n)^F2;
    vector3D T1 = ((-Molecula1.R)*n)^F2;
    Molecula2.AgregueFuerza(F2); Molecula1.AgregueFuerza(F2*(-1));
    Molecula2.AgregueTorque(T2); Molecula1.AgregueTorque(T1);
  }
  else if(Contacto==TRUE){
    ele.cargue(0,0,0);
    Contacto=FALSE;
  }
}

void Colisionador::CalculeTodasLasFuerzas(std::vector<Cuerpo> & Molecula, std::vector<Cuerpo> & Paredes, double dt){
  int i,j,k;// vector3D Gvector;
  
  for(i=0;i<Molecula.size();i++){
    Molecula[i].BorreFuerzaYTorque();
  }

  for(i=0;i<Paredes.size();i++){
    Paredes[i].BorreFuerzaYTorque();
  }

  /*  Gvector.cargue(0,-g,0);
  for(i=0;i<Molecula.size();i++){
    Molecula[i].AgregueFuerza(Gvector*Molecula[i].m);
    }*/
  
  for(i=0;i<Molecula.size();i++){
    for(k=0;k<Paredes.size();k++){
      CalculeFuerzaEntre(Molecula[i],Paredes[k],elePared[i*NumP + k],EstabaEnContactoPared[i*NumP + k],dt);
    }
    for(j=i+1;j<Molecula.size();j++){
      CalculeFuerzaEntre(Molecula[i],Molecula[j],ele[i*N + j],EstabaEnContacto[i*N + j],dt);
    }
  }
}

//---------------- Funciones Globales ----------------

void InicieAnimacion(void){
  std::cout<<"set terminal gif animate"<<std::endl;
  std::cout<<"set output 'dosgranos_2207v5.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange[-10:"<<lx+10<<"]"<<std::endl;
  std::cout<<"set yrange["<<ly/2.0-10<<":"<<ly/2.0+10<<"]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange[0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}

void InicieCuadro(void){
  std::cout<<"plot 0,0 ";
  std::cout<<" , "<<(lx+200)/7<<"*t-100.0,0";
  std::cout<<" , "<<(lx+200)/7<<"*t-100.0,101";
}

void TermineCuadro(void){
  std::cout<<std::endl;
}


int main(){
  
  //Se define una variable de cuerpo llamada Molecula
  std::vector<Cuerpo> Molecula(2);
  std::vector<Cuerpo> Pared(NumP);
  Colisionador Newton;
  Crandom ran64(1);
  double t,tdibujo;
  int tcrear;
  double dt=0.001;
  int i,j;
  double m0 = 1.0, R0 = 1.0;
  double kT = 10.0;
  double V0 = 2.5, deltax;
  double h = 100;

  Newton.Inicie();
  InicieAnimacion();

  //---------- Dos granos -----------
  deltax = 0.02*ran64.r() - 0.01;
  /*    void Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0, double Vz0,
	      double theta0, double omega0,
	      double m0,double R0);*/

  Molecula[0].Inicie(lx/2.0 +5,h/2.0,0,-V0,0,0,3*M_PI/2.0,0,m0,R0); 
  Molecula[1].Inicie(lx/2.0 -5,h/2.0,0,V0,0,0,M_PI/2.0,0,m0,R0); 
  //---------- Paredes ---------------
  double Mpared=100*m0, Rpared = 0.2;
  Cuerpo MolPared;
  for(int i = 0; i < NumP; i++){
    Pared[i].Inicie(-100 + 2*i*Rpared,-Rpared,0,0,0,0,0,0,Mpared,Rpared);
    }

  Cuerpo NuevaMol; 
  double T = (double)N - 1;
    //Se empieza a mover el cuerpo
  for(t=tdibujo=0,tcrear=0;t<T;t+=dt,tdibujo+=dt,tcrear++){

    //----crear una nueva movecula cada cierto tiempo-------
    /*    if(tcrear == 1000){
      deltax = 0.02*ran64.r() - 0.01;
      NuevaMol.Inicie(lx/2.0 +deltax,h,0,0,0,0,0,0,m0,R0);  
      Molecula.push_back(NuevaMol);
      tcrear = 0;
      }*/
    
    if(tdibujo>T/600){ //1200 frames
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
  //std::cout<<Molecula.size()<<'\t'<<Pared.size()<<std::endl;
  return 0;
}
