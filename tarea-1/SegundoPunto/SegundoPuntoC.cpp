#include<iostream>
#include<cmath>
using namespace std;

const double err = 1e-11;

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

double f(double lamb){
	double t,x1,x2,dt=0.0001;
	for(t=0.00001,x1=1,x2=0;t<=1;){
      PasoRK(x1,x2,t,dt,lamb);
    }
	return x1;
}

double bisec(double a, double b){
  double fda,fdx,x0;
  fda = f(a);
  
  while(b-a>err){
    x0 = 0.5*a+0.5*b, fdx = f(x0);
    if(fda*fdx<0){
      b=x0;
    }
    else{
      a=x0,fdx=fda;
    }
  }
  return 0.5*(a+b);
}

void GrafR(double lamb){
  double t,x1,x2,dt=0.001;

  for(t=0.001,x1=1,x2=0;t<1;){
    cout<<t<<" "<<x1<<endl;
    PasoRK(x1,x2,t,dt,lamb);
  }
  cout<<endl;
}

int main(void){
  double a,b,l1,l2,l3,l4,l5;
  
  // a=1.5,b=2.5;
  // l1=bisec(a,b);
  // cout<<"lambda 1 = "<<l1<<endl;
  // GrafR(l1);
  
  // a=5,b=6;
  // l2=bisec(a,b);
  // cout<<"lambda 2 = "<<l2<<endl;
  // GrafR(l2);
  
  // a=8,b=10;
  // l3=bisec(a,b);
  // cout<<"lambda 3 = "<<l3<<endl;
  // GrafR(l3);
  
  // a=11,b=13;
  // l4=bisec(a,b);
  // cout<<"lambda 4 = "<<l4<<endl;
  // GrafR(l4);
  
  a=14,b=16;
  l5=bisec(a,b);
  // cout<<"lambda 5 = "<<l5<<endl;
  GrafR(l5);
  
  return 0;
}
