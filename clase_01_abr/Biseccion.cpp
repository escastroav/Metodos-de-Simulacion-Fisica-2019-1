// Mi primer programa en c++
#include<iostream>
#include<cmath>
using namespace std;

const double err = 1e-10;

double f(double x){
  return sin(x)/x;
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

int main(void){
  double a,b;
  a=1,b=4;
  
 cout<<"f(x)=0 en x="<<bisec(a,b)<<endl;
  return 0;
}
