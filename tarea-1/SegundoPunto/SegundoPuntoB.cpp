#include<iostream>
#include<cmath>
using namespace std;


double f1(double x1, double x2, double t){
  return x2/t;
}

double f2(double x1, double x2, double t,double lamb){
  return -lamb*lamb*x1*t;
}

void PasoRK(double & x1, double & x2, double & t, double dt, double lamb){
  double dx11,dx12,dx13,dx14;
  double dx21,dx22,dx23,dx24;

  dx11=f1(x1,x2,t)*dt;  							dx21=f2(x1,x2,t,lamb)*dt;
  dx12=f1(x1+0.5*dx11,x2+0.5*dx21,t+0.5*dt)*dt;		dx22=f2(x1+0.5*dx11,x2+0.5*dx21,t+0.5*dt,lamb)*dt;
  dx13=f1(x1+0.5*dx12,x2+0.5*dx22,t+0.5*dt)*dt;		dx23=f2(x1+0.5*dx12,x2+0.5*dx22,t+0.5*dt,lamb)*dt;
  dx14=f1(x1+dx13,x2+dx23,t+dt)*dt;  				dx24=f2(x1+dx13,x2+dx23,t+dt,lamb)*dt;

  x1+=(dx11+2*(dx12+dx13)+dx14)/6,x2+=(dx21+2*(dx22+dx23)+dx24)/6,t+=dt;
}

double Flamb(double & x1, double & x2, double & t, double dt, double lamb){
	for(t=0.00001,x1=1,x2=0;t<=1;){
      PasoRK(x1,x2,t,dt,lamb);
    }
	return x1;
}

int main(void){
  double t,x1,x2,dt=0.001,lamb;

  for(lamb = 0.1;lamb<=15.0;lamb+=0.01){
    cout<<lamb<<" "<<Flamb(x1,x2,t,dt,lamb)<<endl;
  }
  
  return 0;
}
