#include <iostream>

int main(){
  std::cout<<"set xlabel \"t\""<<'\n'
	   <<"set ylabel \"Tau\""<<'\n'
	   <<"set xzeroaxis linestyle  1 lt 1 lc -1"<<'\n'
	   <<"plot \"datos1.dat\" u 1:4 w l t \"K=0.1e10\""<<'\n'
	   <<"replot \"datos2.dat\" u 1:4 w l t \"K=0.2e10\""<<'\n'
    	   <<"replot \"datos3.dat\" u 1:4 w l t \"K=0.5e10\""<<'\n'
    	   <<"replot \"datos4.dat\" u 1:4 w l t \"K=1e10\""<<'\n'
    	   <<"replot \"datos5.dat\" u 1:4 w l t \"K=2e10\""<<'\n'
    	   <<"replot \"datos6.dat\" u 1:4 w l t \"K=5e10\""<<'\n'
    	   <<"replot \"datos7.dat\" u 1:4 w l t \"K=10e10\""<<'\n'
	   <<"pause 20"<<std::endl;
}
