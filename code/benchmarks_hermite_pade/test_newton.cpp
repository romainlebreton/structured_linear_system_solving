#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <NTL/ZZX.h>

#include "ZZX_extra.h"
#include "lzz_pX_mosaic_hankel.h"
#include "ZZ_hermite_pade.h"

NTL_CLIENT

// // /*------------------------------------------------------------*/
// /* reads a vecor of polynomials F and an integer vector type  */
// /*------------------------------------------------------------*/
// void read_vecf_type(Vec<ZZX> &vec_F, Vec<long>& type, const string& name){
//   string line_f, line_type;
//   ifstream file;
//   file.open(name);
//   long inp;

//   long nb;
//   if (! getline(file, line_f)){
//     cout << "first line of file should be number of polynomials\n";
//     exit(-1);
//   }
//   istringstream iss_nb{line_f};
//   iss_nb >> nb;
//   vec_F.SetLength(nb);

//   for (long i = 0; i < nb; i++){
//     Vec<ZZ> f;
//     if (! getline(file, line_f)){
//       cout << "not enough polynomials\n";
//       exit(-1);
//     }
//     istringstream iss_f{line_f};
//     while(iss_f >> inp)
//       f.append(ZZ(inp));
//     vec_F[i] = conv<ZZX>(f);
//   }

//   if (! getline(file, line_type)){
//     cout << "last line of file should be type\n";
//     exit(-1);
//   }
//   istringstream iss_type = istringstream(line_type);
//   while (iss_type >> inp)
//     type.append(inp);

//   file.close();
// }


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
int main(int argc, char **argv){

  long n = 3;
  long m = 3;
  long b = 40;

  long prec = n*(m+1)-1;

  Vec<ZZX> f;
  Vec<long> type;
  for(long i = 0; i < n; i++){
    ZZX tmp;
    random(tmp, b, prec);
    f.append(tmp);
    type.append(m);
  }

  hermite_pade_general hp(f, type, prec);
  Vec<ZZX> sol;

  hp.switch_mode(0);
  hp.random_solution(sol);
  cout << sol << endl;

  hp.switch_mode(2);
  hp.random_solution(sol);
  cout << sol << endl;
}
