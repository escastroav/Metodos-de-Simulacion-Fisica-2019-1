//Mi primer programa en C++
#include<iostream>
#include<cmath>
using namespace std;

const double ERR=1e-7;

double f1(double R1,double R2,double r){
  return R2;
}

double f2(double R1,double R2,double lambda,double r){
  return (-R2/r-lambda*lambda*R1);
}

void UnPasoDeRungeKutta4(double & x1,double & x2,double lambda,double & r,double dr){
  double dx11,dx12,dx13,dx14; 
  double dx21,dx22,dx23,dx24; 
  dx11=f1(x1,x2,r)*dr;                              dx21=f2(x1,x2,lambda,r)*dr;
  dx12=f1(x1+0.5*dx11, x2+0.5*dx21, r+0.5*dr)*dr;   dx22=f2(x1+0.5*dx11, x2+0.5*dx21, lambda, r+0.5*dr)*dr;
  dx13=f1(x1+0.5*dx12, x2+0.5*dx22, r+0.5*dr)*dr;   dx23=f2(x1+0.5*dx12, x2+0.5*dx22, lambda, r+0.5*dr)*dr;
  dx14=f1(x1+dx13, x2+dx23, r+dr)*dr;               dx24=f2(x1+dx13, x2+dx23, lambda, r+dr)*dr;

  r+=dr;  x1+=(dx11+2*(dx12+dx13)+dx14)/6;          x2+=(dx21+2*(dx22+dx23)+dx24)/6;
} 

int main(void){
  double r=1,lambda,R1,R2,dr=0.1;

  for(lambda=0.01,R1=1,R2=0;lambda<15;lambda+=0.01 ){
    cout<<lambda<<"\t"<<R1<<endl;
    UnPasoDeRungeKutta4(R1,R2,lambda,r,dr);
  }

  return 0;
}
