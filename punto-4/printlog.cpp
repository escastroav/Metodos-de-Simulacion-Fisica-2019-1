#include <iostream>

int main(){
  std::cout<<"set xlabel \"K\" font \",12\""<<'\n'
	   <<"set ylabel \"{/Symbol t}_{max}\" font \",15\""<<'\n'
	   <<"set logscale xy"<<'\n'
    	   <<"plot \"datos.dat\" u 1:2 w p lc 7 pt 6 ps 2 t \"{/Symbol t}_{max}\""<<'\n'
	   <<"replot \"datos.dat\" u 1:2 w l lc 10 notitle"<<'\n'
	   <<"pause 20"<<std::endl;
}
