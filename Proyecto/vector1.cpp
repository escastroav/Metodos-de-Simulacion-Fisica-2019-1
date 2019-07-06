#include <vector>
#include <iostream>

int main(){
  std::vector <double> Try(5);
  std::cout<<Try.size()<<std::endl;
  int i;
  Try.push_back(i);
  std::cout<<Try.size()<<std::endl;
}
