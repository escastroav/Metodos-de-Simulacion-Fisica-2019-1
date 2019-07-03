#include <iostream>
#include <cmath>
#include "Vector.h"

//------- Define Global Constants -----------------------------------

const int N=3;

const double g=9.8;
const double l=0.12;
const double R=0.015;
const double k=10e10;

const double Zi=0.1786178958448091e0;
const double Lambda=0.2123418310626054*(-1);
const double Xi=0.06626458266981849*(-1);

const double C1=(1-2*Lambda)/2;
const double C2=(1-2*(Xi+Zi));

//------- Define classes ---------------------------------------------

class Body;
class Force;

//------- Body Class -----------------------------------------------

class Body{

private:

  vector3D r, v, F, T;
  double xe;
  double m, R;

public:

  //--- internal functions ----------------------------

  void Start(vector3D r0, vector3D v0, double m, double R, double xi);
  void SetZeroF(void){F.cargue(0,0,0);}
  void SetForce(vector3D Ft){F+=Ft;}
  void SetZeroT(void){T.cargue(0,0,0);}
  void SetTorque(vector3D dT){T+=dT;}
  void SetPosition(double dt, double c);
  void SetVelocity(double dt, double c);
  void Draw(void);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  double Getz(void){return r.z();};
  double GetVx(void){return v.x();};
  double GetVy(void){return v.y();};
  double GetVz(void){return v.z();};
  double GetTx(void){return T.x();};
  double GetTy(void){return T.y();};
  double GetTz(void){return T.z();};
  vector3D GetT(void){return T;};
  vector3D Getr(void){return r;};

  //--- friend functions -------------------------------
  friend class Force;
};

//----- functions Body Class -------------------------

void Body::Start(vector3D r0, vector3D v0, double m0, double R0, double xi){
  r = r0; v = v0; m = m0; R = R0;
  xe = xi;
}

void Body::SetPosition(double dt, double c){
  r += v*(dt*c);
}

void Body::SetVelocity(double dt, double c){
  v += F*(dt*c/m);
}

void Body::Draw(void){
  std::cout<<" , "<<xe<<"+"<<(r.x()-xe)*(1-(R/l))/7<<"*t,"<<r.y()*(1-(R/l))/7<<"*t";
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

//------- Force Class -----------------------------------------------

class Force{
public:
  //--- internal functions ----------------------------
  void SetGravitationAndTensionForce(Body & BodyA);
  void SetForceBTP(Body & BodyA, Body & BodyB); //BTP : Between two Bodies
  void SetTotalForce(Body *Bodies);
  void SetTotalTorque(Body *Bodies);
};

//----- functions Force Class -------------------------

void Force:: SetGravitationAndTensionForce(Body & BodyA){
  vector3D G;

  //----- Calculo tensión de la cuerda ----------------
  double xp = BodyA.Getx() - BodyA.xe;
  vector3D R; R.cargue(xp,BodyA.Gety(),BodyA.Getz());
  double aux = xp*BodyA.GetVy() - BodyA.Gety()*BodyA.GetVx();
  //para calcular la tensión, se tiene en cuenta que el cuerpo esta bajo la acción otra fuerza a parte de la gravedad
  //esa fuerza es la fuerza de Hertz, por lo que es necesario calcularla primero antes de hallar la tensión
  double T = (BodyA.m)*pow(aux,2.0)*pow(l,-3.0) - (BodyA.m)*g*(BodyA.Gety())/l + (BodyA.F*R)/norma(R);

  //----- Calculo Fuerza total sobre el pendulo dada por la tension y por la gravedad-------
  G.cargue(-T*(xp)/l,-(T*(BodyA.Gety())/l+(BodyA.m)*g),0.0);
  BodyA.SetForce(G);
}

void Force::SetForceBTP(Body & BodyA, Body & BodyB){
  vector3D rab = BodyA.r - BodyB.r;
  double s = 2*R - norma(rab);
  vector3D G;

  //------ calculo fuerza de Hertz cuando s > 0 ------------
  if(s > 0){
    G = rab*(k*pow(s,1.5)/norma(rab));
    BodyA.SetForce(G); BodyB.SetForce(G*(-1));
  }
}

void Force::SetTotalForce(Body *Bodies){
  int i,j;
  for(i = 0; i < N; i++){
    Bodies[i].SetZeroF();
  }
  for(i = 0; i < N; i++){
    //calcula primero la fuerza de Hertz entre las esferas, para despues cacular la tension
    for(j = i + 1; j < N; j++){ SetForceBTP(Bodies[i],Bodies[j]); }
    //calcula ahora la fuerza sobre el cuerpo debida a la tensión y a la gravedad
    SetGravitationAndTensionForce(Bodies[i]);
  }
}

void Force::SetTotalTorque(Body *Bodies){
  vector3D Re, R, G;
  int i, j;
  SetTotalForce(Bodies);

  //------- Calculo torque sobe todos los pendulos --------------
  for(i = 0; i < N; i++){
    Bodies[i].SetZeroT();
    Re.cargue(Bodies[i].xe,0,0);
    R = Bodies[i].r - Re;
    G = R^Bodies[i].F;
    Bodies[i].SetTorque(G);
  }

}

//------------- Global Functions -----------------------------------

void StartAnimation(void){
  //std::cout<<"set term gif animate"<<std::endl;
  //std::cout<<"set output \"CunaNewton.gif\""<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange[-0.15:0.15]"<<std::endl;
  std::cout<<"set yrange[-0.15:0.02]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange [0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}

void StartFrame(void){
    std::cout<<"plot 0,0 ";
    std::cout<<" , -0.15 +"<<0.3/7<<"*t,0";

}
void EndFrame(void){
    std::cout<<std::endl;
}

//------------- Main Program Function ------------------------------

int main(){

  Body Bodies[N];
  Force Newton;
  double t, tprint, dt;
  double m,x10,y10,x20,y20,x30,y30,vx0,vy0,theta0;
  double T;
  int i;

  dt = 0.000001; m = 0.1; T = 2*M_PI*sqrt(l/g);

  StartAnimation();

  //---- Definiendo condiciones iniciales de los pendulos ---------
  theta0 =M_PI/12.0;
  x10 = 0.0; y10 = -l;
  x20 = 2*R; y20 = -l;
  x30 = -2*R - l*sin(theta0); y30 = -l*cos(theta0);
  vx0 = 0.0; vy0 = 0.0;

  vector3D X10, V0;
  vector3D X20;
  vector3D X30;

  X10.cargue(x10,y10,0); V0.cargue(vx0,vy0,0);
  X20.cargue(x20,y20,0);
  X30.cargue(x30,y30,0);

  Bodies[0].Start(X10,V0,m,R,0.0);
  Bodies[1].Start(X20,V0,m,R,2*R);
  Bodies[2].Start(X30,V0,m,R,-2*R);

  for(t=tprint=0; t<5*T; t+=dt, tprint+=dt){
    if(tprint > T/100){
      StartFrame();
      for(i=0;i<N;i++){Bodies[i].Draw();}
      EndFrame();
      tprint=0;
      }

    //------- exporta datos torque sobre pendulo del medio ---------------------------------
    /*if(norma(Bodies[0].Getr()-Bodies[1].Getr())<2*R || norma(Bodies[0].Getr()-Bodies[2].Getr())<2*R){
      Newton.SetTotalTorque(Bodies);
      std::cout<<t<<'\t'<<Bodies[0].GetTz()<<std::endl;
      }*/

    for(i=0;i<N;i++) Bodies[i].SetPosition(dt,Zi);
    Newton.SetTotalForce(Bodies);
    for(i=0;i<N;i++) Bodies[i].SetVelocity(dt,C1);
    for(i=0;i<N;i++) Bodies[i].SetPosition(dt,Xi);
    Newton.SetTotalForce(Bodies);
    for(i=0;i<N;i++) Bodies[i].SetVelocity(dt,Lambda);
    for(i=0;i<N;i++) Bodies[i].SetPosition(dt,C2);
    Newton.SetTotalForce(Bodies);
    for(i=0;i<N;i++) Bodies[i].SetVelocity(dt,Lambda);
    for(i=0;i<N;i++) Bodies[i].SetPosition(dt,Xi);
    Newton.SetTotalForce(Bodies);
    for(i=0;i<N;i++) Bodies[i].SetVelocity(dt,C1);
    for(i=0;i<N;i++) Bodies[i].SetPosition(dt,Zi);
  }
}
