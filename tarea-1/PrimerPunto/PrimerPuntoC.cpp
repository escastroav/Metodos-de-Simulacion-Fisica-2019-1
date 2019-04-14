#include<iostream>
#include<cmath>
using namespace std;

const double beta2 = 0.08;

double f1(double x1, double x2, double t, double beta){
  return -beta*x1*x2;
}

double f2(double x1, double x2, double t,double beta){
  return beta*x1*x2-beta2*x2;
}

void PasoRK(double & x1, double & x2, double & t, double dt,double beta){
  double dx11,dx12,dx13,dx14;
  double dx21,dx22,dx23,dx24;

  dx11=f1(x1,x2,t, beta)*dt;  							dx21=f2(x1,x2,t, beta)*dt;
  dx12=f1(x1+0.5*dx11,x2+0.5*dx21,t+0.5*dt, beta)*dt;	dx22=f2(x1+0.5*dx11,x2+0.5*dx21,t+0.5*dt, beta)*dt;
  dx13=f1(x1+0.5*dx12,x2+0.5*dx22,t+0.5*dt,beta)*dt;  	dx23=f2(x1+0.5*dx12,x2+0.5*dx22,t+0.5*dt,beta)*dt;
  dx14=f1(x1+dx13,x2+dx23,t+dt,beta)*dt;  				dx24=f2(x1+dx13,x2+dx23,t+dt,beta)*dt;

  x1+=(dx11+2*(dx12+dx13)+dx14)/6,x2+=(dx21+2*(dx22+dx23)+dx24)/6,t+=dt;
}

int main(void){
  double t,x1,x2,x3,dt=0.01,beta;

  for(beta=0.0;beta<0.6;beta+=0.001){
    for(t=0,x1=0.999,x2=0.001,x3=0.;t<10000;){
      PasoRK(x1,x2,t,dt,beta);
      if(t-9999.5>0){
	cout<<beta/beta2<<" "<<x1<<endl;
	break;
      }
    }
  }
  
  
  
  return 0;
}
