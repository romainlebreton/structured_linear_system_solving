#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;

/*************************************************
* Generates test cases where there are t entries *
* of n with rank defect 1                        *
*************************************************/
int main(){
  for (long nbits = 1; nbits < 30; nbits++){
    for (long n = 1; n < 60; n++){
      for (long t = 1; t <= 5; t++){
        ostringstream oss;
        oss << "test" << nbits << "_" << t << "x" << n << ".data";
        string filename = oss.str();
        cout << "creating: " << filename << endl;
        ofstream ofs{filename};
        ofs << nbits << endl;
        ofs << 1 << endl;
        for (long i = 0; i < n; i++){
	  ofs << t-1 << " ";
      	}
      }
    }
  }
}
