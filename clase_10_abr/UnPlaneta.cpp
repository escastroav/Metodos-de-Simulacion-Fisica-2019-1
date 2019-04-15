#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double GM = 1.0;

//Definicion de un cuerpo como molecula
class Cuerpo
{
  
private:
  double
  x,
    y,
    Vx,
    Vy,
    Fx,
    Fy,
    m,
    R;
  
public:
  
  //Ordenes
  void Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0);
  void CalculeFuerza(void);
  void Muevase(double dt);

  //Impresion
  void Dibujese(void);

  //Lectura de datos
  double Getx(void){return x;}; //Funcion Inline, para evitar que el compilador se muevea  dos puntos de la memoria
  double Gety(void){return y;};
}; //Ojo! este punto y coma es obligatorio!!

void Cuerpo::Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0)
{
  x = x0;
  y = y0;
  Vx = Vx0;
  Vy = Vy0;
  m = m0;
  R = R0;  
}

void Cuerpo::CalculeFuerza(void)
{
  double aux = GM*m*pow(x*x+y*y,-1.5);
  Fx = -aux*x;
  Fy = -aux*y;
}

void Cuerpo::Muevase(double dt)
{
  x += Vx*dt;
  y += Vy*dt;

  Vx += Fx/m * dt;
  Vy += Fy/m * dt;
}

void Cuerpo::Dibujese(void)
{
  cout << " , " << x << "+" << R << "*cos(t), " << y << "+" << R << "*sin(t)";
}

//Funciones de impresion

void InicieAnimacion(void)
{
  //cout << "set terminal gif animate" << endl;
  //cout << "set output 'balon.gif'" << endl;
  cout << "unset key" << endl;
  cout << "set xrange [-110:110]" << endl;
  cout << "set yrange [-110:110]" << endl;
  cout << "set size ratio -1" << endl;
  cout << "set parametric" << endl;
  cout << "set trange [0:7]" << endl;
  cout << "set isosamples 12" << endl;
}

void InicieCuadro(void)
{
  cout << "plot 0,0 ";
}
void TermineCuadro(void)
{
  cout << endl;
}

int main()
{
  Cuerpo Planeta;
  double t,tdibujo,dt = 1.0;
  double r=50,omega,T,V;

  InicieAnimacion();

  omega=sqrt(GM/(r*r*r));
  T=2*M_PI/omega;
  V=omega*r;

  //Dar orden de imponer las condiciones iniciales (x,y,Vx,Vy,m,R)
  Planeta.Inicie(r, 0, -V*0.1, V, 1.0, 5);

  //Evolucion del balon
  for(t=tdibujo=0; t<1.1*T; t+=dt, tdibujo+=dt)
    {
      //cout << Planeta.Getx() << "\t" << Planeta.Gety() << endl;
      if(tdibujo>T/1000)
	{
	  InicieCuadro();
	  Planeta.Dibujese();
	  TermineCuadro();
	  tdibujo=0;
	}
      Planeta.CalculeFuerza();
      Planeta.Muevase(dt);
    }
  
  return 0;
}
