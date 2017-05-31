#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;

/*************************************************
* Generates test cases where there are t entries *
* of n with rank defect 1                        *
*************************************************/
int main(){
  long t = 5;
  long nbits = 10;
  for (long q = 1; q <= 20; q++){
    long n = q * 50;
    ostringstream oss;
    oss << "test" << nbits << "_" << n << "x" << t << ".data";
    string filename = oss.str();
    cout << "creating: " << filename << endl;
    ofstream ofs{filename};
    ofs << nbits << endl;
    ofs << 1 << endl;
    for (long i = 0; i < t; i++){
       ofs << n-1 << " ";
    }
  }
}

