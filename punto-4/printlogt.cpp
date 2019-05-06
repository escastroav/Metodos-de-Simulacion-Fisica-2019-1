#include <iostream>

int main(){
  std::cout<<"set xlabel \"K\" font \",12\""<<'\n'
	   <<"set ylabel \"t_{max}\" font \",12\""<<'\n'
	   <<"set logscale xy"<<'\n'
    	   <<"plot \"datos.dat\" u 1:3 w p lc 9 pt 6 ps 2 t \"t_{max}\""<<'\n'
	   <<"replot \"datos.dat\" u 1:3 w l lc 6 notitle"<<'\n'
	   <<"pause 20"<<std::endl;
}
