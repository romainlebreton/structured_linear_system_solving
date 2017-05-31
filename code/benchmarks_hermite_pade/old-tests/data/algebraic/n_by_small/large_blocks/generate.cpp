#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;

/*************************************************
* Generates test cases where there are n entries *
* of n with rank defect 1                        *
*************************************************/
int main(){
  long nbits = 1000;
  for (long dy = 5; dy <= 10; dy += 5)
    for (long dx = 10; dx <= 490; dx += 100){
      ostringstream oss;
      oss << "test" << nbits << "_" << dx << "x" << dy << ".data";
      string filename = oss.str();
      cout << "creating: " << filename << endl;
      ofstream ofs{filename};
      ofs << nbits << endl;
      ofs << dx << endl;
      ofs << dy << endl;
    }
}
