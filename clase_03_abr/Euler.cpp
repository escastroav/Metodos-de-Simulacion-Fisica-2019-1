#include<iostream>
#include<cmath>
using namespace std;

const double err = 1e-10;

double f(double x, double t){
  return x;
}

void PasoEuler(double & x, double & t, double dt){
  double dx;

  dx=f(x,t)*dt;
  x+=dx,t+=dt;
}

int main(void){
  double t,x,dt=0.01;

  for(t=0,x=1;t<2;){
    cout<<t<<" "<<x<<endl;
    PasoEuler(x,t,dt);
  }
  
  return 0;
}
