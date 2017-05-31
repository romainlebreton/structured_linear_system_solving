#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;

/*************************************************
* Generates test cases where there are n entries *
* of n with rank defect 1                        *
*************************************************/
int main(){
  long nbits = 10;
  for (long t = 1; t < 20; t++){
    long n = t * 10;
    ostringstream oss;
    oss << "test" << nbits << "_" << n << ".data";
    string filename = oss.str();
    cout << "creating: " << filename << endl;
    ofstream ofs{filename};
    ofs << nbits << endl;
    ofs << 1 << endl;
    for (long i = 0; i < n; i++){
      ofs << n-1 << " ";
    }
  }
}
