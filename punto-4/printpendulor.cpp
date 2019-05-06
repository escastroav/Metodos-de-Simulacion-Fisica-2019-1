#include<iostream>

int main(){
  std::cout<<"set xlabel \"(t-t0)K^{-b}\" font \",12\""<<'\n'
	   <<"set ylabel \"{/Symbol t}K^{-a}\" font \",12\""<<'\n'
	   <<"plot \"datos1r.dat\" w l t \"K=0.1e10\""<<'\n'
	   <<"replot \"datos2r.dat\" w l t \"K=0.2e10\""<<'\n'
	   <<"replot \"datos3r.dat\" w l t \"K=0.5e10\""<<'\n'
	   <<"replot \"datos4r.dat\" w l t \"K=1e10\""<<'\n'
	   <<"replot \"datos5r.dat\" w l t \"K=2e10\""<<'\n'
	   <<"replot \"datos6r.dat\" w l t \"K=5e10\""<<'\n'
	   <<"replot \"datos7r.dat\" w l t \"K=10e10\""<<'\n'
	   <<"pause 20"<<std::endl;

}
