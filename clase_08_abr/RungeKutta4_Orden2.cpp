//Mi primer programa en C++
#include<iostream>
#include<cmath>
using namespace std;

const double ERR=1e-7;

const double omega=1.0;

double f1(double x1,double x2,double t){
  return x2;
}

double f2(double x1,double x2,double t){
  return -omega*omega*x1;
}

void UnPasoDeRungeKutta4(double & x1,double & x2,double & t,double dt){
  double dx11,dx12,dx13,dx14; 
  double dx21,dx22,dx23,dx24; 
  dx11=f1(x1,x2,t)*dt;                              dx21=f2(x1,x2,t)*dt;
  dx12=f1(x1+0.5*dx11, x2+0.5*dx21, t+0.5*dt)*dt;   dx22=f2(x1+0.5*dx11, x2+0.5*dx21, t+0.5*dt)*dt;
  dx13=f1(x1+0.5*dx12, x2+0.5*dx22, t+0.5*dt)*dt;   dx23=f2(x1+0.5*dx12, x2+0.5*dx22, t+0.5*dt)*dt;
  dx14=f1(x1+dx13, x2+dx23, t+dt)*dt;               dx24=f2(x1+dx13, x2+dx23, t+dt)*dt;

  t+=dt;  x1+=(dx11+2*(dx12+dx13)+dx14)/6;          x2+=(dx21+2*(dx22+dx23)+dx24)/6;
} 

int main(void){
  double t,x1,x2,dt=0.1;

  for(t=0,x1=0,x2=1;t<10; ){
    cout<<t<<" "<<x1<<endl;
    UnPasoDeRungeKutta4(x1,x2,t,dt);
  }

  return 0;
}
