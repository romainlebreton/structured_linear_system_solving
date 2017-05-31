#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;

/*************************************************
* Generates test cases where there are n entries *
* of n with rank defect 1                        *
*************************************************/
int main(){
  for (long nbits = 100; nbits <= 2000; nbits += 100){
    for (long m = 1; m <= 10; m++){
      long n = m*10;
      ostringstream oss;
      oss << "test" << nbits << "_" << n << ".data";
      string filename = oss.str();
      cout << "creating: " << filename << endl;
      ofstream ofs{filename};
      ofs << nbits << endl;
      ofs << n << endl;
      ofs << n << endl;
    }
  }
}
