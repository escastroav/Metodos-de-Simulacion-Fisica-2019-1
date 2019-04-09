#include<iostream>
#include<cmath>
using namespace std;

const double err = 1e-10;

double f(int n,double x, double t){
  return cos(n*t-x*sin(t));
}

double Simpson(double a, double b, int N, int n, double x){
  double h=(b-a)/(2*N),t,suma=0;
  
  suma+=f(n,x,a)+f(n,x,a+h*2*N);
  
  for(int i=1;i<=N;i++){
    t=a+h*2*i;
    suma+=2*f(n,x,t);
  }
  
  for(int i=1;i<=N;i++){
    t=a+h*(2*i-1);
    suma+=4*f(n,x,t);
  }
  
  return suma*h/3.0;
}

double Bessel(int n, double x){
  int NN = 500;
  return 1/M_PI*Simpson(0,M_PI,NN,n,x);
}

int main(void){
  double x;
  int n=0;

  for(x=0;x<10;x+=0.01)
    cout<<x<<" "<<Bessel(n,x)<<endl;
  
  return 0;
}
