#include<iostream>
#include<sstream>
#include<vector>
#include<fstream>
#include<algorithm>
using namespace std;

struct result{
  double total = 0;
  double mul_total = 0, m_right = 0, m_left = 0;
  double x_right = 0, x_left = 0, y_right = 0, y_left = 0;
  double mul_lower = 0;
  double cauchy_left = 0, cauchy_right = 0;
  double recon = 0;
  long rank = 0, nbits = 0;
  
  static const int N_STR = 15;

  string str[N_STR] = 
    {"rank", "nbits","everything", "mul time", "mul M right", 
     "mul X right", "mul Y right", "mul M left", "mul X left",
     "mul Y left", "cauchy left", "cauchy right",
     "total reconstruction time", "for lower", "Took"};

  int search(string s){
    for (int i = 0; i < N_STR; i++){
      int index = s.find(str[i]);
      if (index != string::npos) return i;
    }
    return -1;
  }

  double get_value(string s){
    istringstream iss{s};
    double val = 0;
    while (true){
      iss >> val;
      if (!iss){
        if (iss.eof()) break;
	else{
	  iss.clear();
	  iss.ignore();
	}
      }
    }
    return val;
  }

  result(ifstream &ifs){
    string line;
    while (getline(ifs, line)){
      int t = search(line);
      if (t == -1) continue;
      double val = get_value(line);
      //cout << "val: " << val << endl;
      switch(t){
        case  0: rank = (long)val;
                break;
        case  1: nbits = (long)val;
                break;
        case  2: if (total == 0) total = val;
                break;
        case  3: mul_total = val;
                break;
        case  4: m_right = val; break;
	case  5: x_right = val; break;
	case  6: y_right = val; break;
	case  7: m_left = val; break;
	case  8: x_left = val; break;
	case  9: y_left = val; break;
	case 10: cauchy_left = val; break;
	case 11: cauchy_right = val; break;
	case 12: recon = val; break;
	case 13: mul_lower = val; break;
	case 14: if (total == 0) total = val; break;
      }
    }
  }

  void print(){
    cout << rank << " " << nbits << " " << total << " " << mul_total
    << " " << m_right << " " << x_right << " " << y_right << " "
    << m_left << " " << x_left << " " << y_left << " " << cauchy_left
    << " " << cauchy_right << " " << recon << " " << mul_lower << endl;
  }
};

bool compare (const result &i, const result &j){
  return i.nbits < j.nbits;
}

int main(int argc, char* argv[]){
  vector<result> rs;  
  for(long i = 1; i < argc; i++){
    ifstream ifs{argv[i]};
    result r{ifs};
    rs.emplace_back(r);
  }
  sort(rs.begin(), rs.end(), compare);
  for (auto &i: rs) i.print();
}

















