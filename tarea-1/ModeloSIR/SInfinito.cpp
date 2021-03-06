//Mi primer programa en C++
#include<iostream>
#include<cmath>
using namespace std;

const double gma=0.08;

double fs(double s,double i, double r, double btagma, double t){
  return -btagma*gma*s*i;
}

double fi(double s,double i, double r, double btagma, double t){
  return (btagma*gma*s*i - gma*i);
}

double fr(double s,double i, double r,double t){
  return gma*i;
}

void UnPasoDeRungeKutta4(double & s,double & i,double & r, double btagma, double & t,double dt){
  double ds1,ds2,ds3,ds4; 
  double di1,di2,di3,di4;
  double dr1,dr2,dr3,dr4;
  
  ds1=fs(s,i,r,btagma,t)*dt;
  di1=fi(s,i,r,btagma,t)*dt;
  dr1=fr(s,i,r,t)*dt;
  
  ds2=fs(s+0.5*ds1, i+0.5*di1, r+0.5*dr1, btagma, t+0.5*dt)*dt;
  di2=fi(s+0.5*ds1, i+0.5*di1, r+0.5*dr1, btagma, t+0.5*dt)*dt;
  dr2=fr(s+0.5*ds1, i+0.5*di1, r+0.5*dr1, t+0.5*dt)*dt;

  ds3=fs(s+0.5*ds2, i+0.5*di2, r+0.5*dr2, btagma, t+0.5*dt)*dt;
  di3=fi(s+0.5*ds2, i+0.5*di2, r+0.5*dr2, btagma, t+0.5*dt)*dt;
  dr3=fr(s+0.5*ds2, i+0.5*di2, r+0.5*dr2,  t+0.5*dt)*dt;
 
  ds4=fs(s+ds3, i+di3, r+dr3, btagma, t+dt)*dt;
  di4=fi(s+ds3, i+di3, r+dr3, btagma, t+dt)*dt;
  dr4=fr(s+ds3, i+di3, r+dr3, t+dt)*dt;

  t+=dt;
  s+=(ds1+2*(ds2+ds3)+ds4)/6;
  i+=(di1+2*(di2+di3)+di4)/6;
  r+=(dr1+2*(dr2+dr3)+dr4)/6;
} 

int main(void){
  double t,s,i,r,beta,dt=0.1;

  for(t=1000,s=0.999,i=0.001,r=0,beta=0.08;beta<2.0;beta+=0.002){
    cout<<beta/gma<<"\t"<<s<<endl;
    UnPasoDeRungeKutta4(s,i,r,beta/gma,t,dt);
  }

  return 0;
}
