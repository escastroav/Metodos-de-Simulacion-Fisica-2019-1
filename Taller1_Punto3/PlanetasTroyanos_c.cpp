#include<iostream>
#include<cmath>
#include<fstream>
#include "Vector.h"

const int N = 3;

const double G=1.0;

//Constantes para PEFRL
const double Zhi = 0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Chi = -0.06626458266981849;
const double UnoMenosDosLambdaDivididoDos = (1.0-2*Lambda)/2.0;
const double UnoMenosDosChiMasZhi = 1-2*(Zhi+Chi);



//-------------Declaración de las clases------------------
class Cuerpo;
class Colisionador;

//------------Clase Cuerpo-------------------------------

class Cuerpo{
  //Los datos en cuerpo son privados, las ordenes son publicas
private:     
  vector3D r,V,F; double m,R;

public:
  void Inicie(double x0,double y0, double z0, double Vx0,double Vy0, double Vz0, double m0, double R0);
  void BorreFuerza(void){F.cargue(0,0,0);};//Inline
  void AgregueFuerza(vector3D dF){F+=dF;};//Inline, le agregamos a la fuerza un vector que pasamos

  //Usando el método de Forest-Ruth, necesito devolverme y calcular la posición y la velocidad usando la velocidad y la fuerza/m respectivamente
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);

  //Se hacen funciones separadas macro, i.e, que vaya haciendo la funcion mientras se compila el programa, para sacar los datos privados y solo imprimirlos y no modificarlos
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();}; //Inline

  //Le doy permiso a la clase colisionador de que acceda a alos datos de Cuerpo
  friend class Colisionador;
};

//Para que el compilador sepa que una funcion es de la clase cuerpo, se le escribe lo del principio, y puede usar las variables de la clase cuerpo

//Se definen las condiciones iniciales del cuerpo
void Cuerpo::Inicie(double x0,double y0, double z0, double Vx0,double Vy0,double Vz0, double m0,double R0){
  r.cargue(x0,y0,z0);//Se cargan los vectores con la posición y velocidad iniciales
  V.cargue(Vx0,Vy0,Vz0);
  m=m0; R=R0;//Y se definen masa y radio iniciales.
}


void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(dt*Coeficiente);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(dt*Coeficiente/m);//Dividimos por la masa para usar la aceleración
}


//Para dibujar el balón
void Cuerpo::Dibujese(void){//Es la instrucción para gnuplot de dibujar puntos en las esquinas del círculo
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}


//-----------------Clase Colisionador----------------
class Colisionador{
private:
public:
  void CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2);
  void CalculeTodasLasFuerzas(Cuerpo * Planetas);//* Planeta porque quiero que me manden una matríz de planetas.
};

//Calculamos la fuerza siguiendo \vec(F) = Gm1m2/r^2 \vec(r_unitario)
void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2){
  vector3D dr = Planeta2.r-Planeta1.r;
  double aux = G*Planeta1.m*Planeta2.m*pow(norma2(dr),-1.5);
  vector3D F1 = dr*aux;
  Planeta1.AgregueFuerza(F1); Planeta2.AgregueFuerza(F1*(-1));
}



void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Planetas){
  int i,j;
  for(i=0; i<N; i++) Planetas[i].BorreFuerza();
  for(i=0; i<N; i++)
    for(j=i+1; j<N;j++)
      CalculeFuerzaEntre(Planetas[i], Planetas[j]);
}

//--------------Funciones Globales ------------------
void InicieAnimacion(void){
  //std::cout<<"set terminal gif animate"<<std::endl;
  //std::cout<<"set output 'pelicula.gif'"<<std::endl;//El nombre del gif
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange[-110:110]"<<std::endl;
  std::cout<<"set yrange[-110:110]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl;//Una curva paramétrica
  std::cout<<"set trange[0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;//La densidad de lineas
}
//Para animar se dibuja un cuadro y luego se dibuja otro encima. Para ello debe iniciar el cuadro, luego retornar y luego plotear sobre el cuadro anterior

void InicieCuadro(){
  std::cout<<"plot 0,0 ";//Inicia con un punto en 0,0

}
void TermineCuadro(){
  std::cout<<std::endl;//Termina la linea escrita en el paso anterior.
}

int main(){
  Cuerpo Planeta[N];//Tenemos muchos planetas, entonces defino una matríz de N planetas
  Colisionador Newton; //El que calcula la fuerza.
  double t, dt=1.0, tdibujo;
  double r = 1000, omega, T;//Se definen distancia entre los planetas, frecuencia de orbita, periodo
  double m0=1047, m1=1, m3 = 0.005;
  double x0, x1, V0, V1, M;
  double x_rot, y_rot;
  int i;
  double c60=cos(M_PI/3), s60=sin(M_PI/3);

  //  InicieAnimacion(); //Primero iniciamos gnuplot

  M = m0+m1+m3; omega = sqrt((G*M)/(r*r*r)); T = 2*M_PI/omega;
  x1 = m1*r/(m0+m1+m3*c60); x0=x1-r; V0=omega*x0; V1=omega*x1;

  //Se dan condiciones iniciales a los planetas
  //----------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
  Planeta[0].Inicie(x0,0,0, 0, V0, 0, m0,10);
  Planeta[1].Inicie(x1,0,0, 0, V1, 0, m1,5);
  Planeta[2].Inicie(x1*c60,x1*s60,0,0,V1,0,m3,1);
  
  //Se empieza a mover el cuerpo
  for(t=tdibujo=0;t<20*T;t+=dt,tdibujo+=dt){
    //Se imprimen los datos a partir de las funciones Get
    //    std::cout<<Planeta[1].Getx()<<"\t"<<Planeta[1].Gety()<<std::endl;//Coord de Júpiter

    //    x_rot=cos(omega*t)*Planeta[1].Getx()+sin(omega*t)*Planeta[1].Gety();
    //    y_rot=-sin(omega*t)*Planeta[1].Getx()+cos(omega*t)*Planeta[1].Gety();
    //    std::cout<<x_rot<<" "<<y_rot<<std::endl;//Coord del sistema que rota con júpiter
    std::cout<<Planeta[3].Getx()<<"\t"<<Planeta[3].Gety()<<std::endl;//Coord de Júpiter

    
    //Para imprimir un número definido de frames se revisa que tdibujo sea un múltiplo de algo
    /*if(tdibujo>T/1000){//el denominador es el número de frames
      InicieCuadro();
      for(i=0; i<N; i++) Planeta[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
      }*/
    //Mover según PEFRL Orden 4, como hay varios planetas debemos moverlos todos
    for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,Zhi);
    
    Newton.CalculeTodasLasFuerzas(Planeta);//Calculamos la aceleración en ese sitio y con ésta calculamos V
    for(i=0; i<N; i++) Planeta[i].Mueva_V(dt, UnoMenosDosLambdaDivididoDos);//Con ésta V calculamos r
    
    for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,Chi);
    
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(i=0; i<N; i++) Planeta[i].Mueva_V(dt,Lambda);
    
    for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,UnoMenosDosChiMasZhi);
    
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(i=0; i<N; i++) Planeta[i].Mueva_V(dt,Lambda);
    
    for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,Chi);

    Newton.CalculeTodasLasFuerzas(Planeta);
    for(i=0; i<N; i++) Planeta[i].Mueva_V(dt, UnoMenosDosLambdaDivididoDos);
    
    for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,Zhi);

  }
  return 0;
}
