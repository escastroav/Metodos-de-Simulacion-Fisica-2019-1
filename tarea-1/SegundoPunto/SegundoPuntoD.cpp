#include<iostream>
#include<cmath>
using namespace std;

const double err = 1e-11;

//---------- Funcion a integrar -------

double f(int n,double x, double t){
  return cos(n*t-x*sin(t));
}

//---------- Metodo de Simpson --------

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

//---------- Funcion de Bessel --------

double Bessel(int n, double x){
  int NN = 100000;
  return 1/M_PI*Simpson(0,M_PI,NN,n,x);
}

//---------- Metodo de Biseccion ------

double bisec(double a, double b, int n){
  double fda,fdx,x0;
  fda = Bessel(n,a);
  
  while(b-a>err){
    x0 = 0.5*a+0.5*b, fdx = Bessel(n,x0);
    if(fda*fdx<0){
      b=x0;
    }
    else{
      a=x0,fdx=fda;
    }
  }
  return 0.5*(a+b);
}

//----------- main --------------------

int main(void){
  double x;
  int n=0;

  /*for(x=0;x<15;x+=0.01)
    cout<<x<<" "<<Bessel(n,x)<<endl;*/

  cout<<"lambda 1 = "<<bisec(2.,3.,n)<<endl;
  cout<<"lambda 2 = "<<bisec(5.0,6.0,n)<<endl;
  cout<<"lambda 3 = "<<bisec(8.0,9.0,n)<<endl;
  cout<<"lambda 4 = "<<bisec(11.0,12.0,n)<<endl;
  cout<<"lambda 5 = "<<bisec(14.0,15.0,n)<<endl;
  
  
  return 0;
}

