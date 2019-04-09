#include<iostream>
#include<cmath>
using namespace std;


/*Resuelve la integral de una funcion por simpson
entre a y b con N puntos
 */

const double err = 1e-10;

double f(double x){
  return sin(x);
}

double Simpson(double a, double b, int N){
  double h=(b-a)/(2*N),x,suma=0;
  suma+=f(a)+f(a+h*2*N);
  
  for(int i=1;i<=N;i++){
    x=a+h*2*i;
    suma+=2*f(x);
  }
  
  for(int i=1;i<=N;i++){
    x=a+h*(2*i-1);
      suma+=4*f(x);
  }
  
  return suma*h/3.0;
}

int main(void){
  double a,b;
  int N=5;
  a=0,b=M_PI;
  
  // cout<<"f(x)=0 en x="<<bisec(a,b)<<endl;
  cout<<"la integral da: "<<Simpson(a,b,N)<<" con "<<N<<" puntos"<<endl;
  return 0;
}
