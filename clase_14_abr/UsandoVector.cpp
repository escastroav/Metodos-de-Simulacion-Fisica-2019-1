#include <iostream>
#include <cmath>
#include "Vector.h"

using namespace std;

int main(void)
{
  vector3D a,b,c;
  double p;

  a.cargue(1,2,0);
  b.cargue(1,1,0);

  p = a*b;
  
  cout << p << endl;
  //c.show();
  return 0;
}
