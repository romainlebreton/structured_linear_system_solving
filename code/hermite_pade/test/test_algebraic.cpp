#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <NTL/ZZX.h>

#include "magma_output.h"
#include "lzz_pX_mosaic_hankel.h"
#include "ZZX_extra.h"
#include "ZZXY.h"
#include "ZZ_hermite_pade.h"
#include "ZZ_pX_extra.h"
#include "ZZ_pXY.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates a polynomial F                                     */
/*  deg(F,x)=deg(F,y)=d, ht(F)=b                              */
/* returns the solution with d^2 terms                        */
/*------------------------------------------------------------*/
void random_F_and_solution(ZZXY & F, ZZX & g, ZZ& den_g, long d, long b){
  random(F, b, d, d);
  ZZX F0;
  random(F0, b, F.degY()-1);
  ZZX factor;
  ZZ root{22};
  SetCoeff(factor, 1, 1);
  SetCoeff(factor, 0, -root);
  F0 = F0 * factor;
  for (long i = 0; i <= F.degY(); i++)
    F.coeffX[i] = (F.coeffX[i] << 1) + coeff(F0, i);
  F.series_solution(g, den_g, root, 4*d*d);
}

/*------------------------------------------------------------*/
/* creates a polynomial F                                     */
/*  deg(F,x)=deg(F,y)=d, ht(F)=b                              */
/* returns the solution with d^2 terms                        */
/*------------------------------------------------------------*/
void random_F_and_solution_mod_p(ZZXY & F, ZZX & g, long dx, long dy, long b){
  zz_p::FFTInit(0);
  ZZ p{zz_p::modulus()};
  // cout << "original prime: " << p << endl;
  while (p < (ZZ{1} << (2*b)))
    p = p*p;
  p = p*p;
  // cout << "original prime after power: " << p << endl;
  ZZ_p::init(p);

  random(F, b, dx, dy);

  ZZX F0;
  random(F0, b, F.degY()-1);
  ZZX factor;
  ZZ root{22};
  SetCoeff(factor, 1, 1);
  SetCoeff(factor, 0, -root);
  F0 = F0 * factor;

  for (long i = 0; i <= F.degY(); i++)
    F.coeffX[i] = (F.coeffX[i] << 1) + coeff(F0, i);

  ZZ_pXY Fp;
  conv(Fp, F);
  
  ZZ_pX g_p;

  Fp.series_solution(g_p, to_ZZ_p(root), 2*(dx+3)*(dy+3));
  conv(g, g_p);
}

/*------------------------------------------------------------*/
/* reads a polynomial f and an integer vector                 */
/*------------------------------------------------------------*/
void read_f_type(ZZX &F, Vec<long>& type, const string& name){
  string line_f, line_type;
  ifstream file;
  file.open(name);
  long inp;
  Vec<ZZ> f;

  if (! getline(file, line_f)){
    cout << "first line of file should be coefficients of f\n";
    exit(-1);
  }
  istringstream iss_f{line_f};
  while(iss_f >> inp)
    f.append(ZZ(inp));
  F = conv<ZZX>(f);
    
  if (! getline(file, line_type)){
    cout << "first line of file should be type\n";
    exit(-1);
  }
  istringstream iss_type = istringstream(line_type);
  while (iss_type >> inp)
    type.append(inp);

  file.close();
}

/*------------------------------------------------------------*/
/* single test. mode_lifting is as in ZZ_hermite_pade.h       */
/* mode_stop = 0 for a witness solution                       */
/* otherwise, waits until we see twice the same result        */
/*------------------------------------------------------------*/
void do_test(long b, long dx, long dy, long mode_lifting, long mode_stop){

  long p2;
  zz_p::FFTInit(1);
  p2 = zz_p::modulus();

  Vec<long> type;
  ZZXY F;
  ZZX f;
  random_F_and_solution_mod_p(F, f, dx, dy, b);
  type.SetLength(dy+2);
  for (long i = 0; i < dy+2; i++)
    type[i] = dx+2;
  
  hermite_pade_algebraic hp(f, type, deg(f)+1, 0);
  if (hp.NumCols() - hp.Rank() > 1){
    Vec<long> new_type = hp.find_new_type();
    hp.init(new_type);
    type = new_type;
  }

  hermite_pade_algebraic hp2(f, type, hp.Rank()+1, 0);
  long sigma = deg(f)+1;
  cout << "rk=" << hp.Rank() << " rk2=" << hp2.Rank() << endl;
  if (hp.Rank() == hp2.Rank())
    sigma = hp2.Rank();
  
  Vec<long> witness;
  zz_p::init(p2);
  for (long i = 0; i < type.length(); i++)
    for (long j = 0; j < type[i]+1; j++)
      if (i <= F.degY())
	witness.append(conv<long>( conv<zz_p>(coeff(F.coeffX[i], j)) / conv<zz_p>(coeff(F.coeffX[0], deg(F.coeffX[0]))) ));
      else
	witness.append(0);
  
  if (mode_stop == 0)
    hp = hermite_pade_algebraic(f, type, sigma, witness, 0);
  else 
    hp = hermite_pade_algebraic(f, type, sigma, 0);
  
  hp.switch_mode(mode_lifting);

  cout << "dx=" << dx << " dy=" << dy << endl;
  cout << "type " << type << endl;
  cout << "dimensions:= " << hp.NumRows() << " " << hp.NumCols() << endl;
  cout << "nullity " << hp.NumCols() - hp.Rank() << endl;
  
  Vec<ZZX> sol;
  hp.random_solution(sol);
  cout << "first coefficient of sol: " << coeff(sol[0], 0) << endl;
 }

/*------------------------------------------------------------*/
/* usage:                                                     */
/*------------------------------------------------------------*/
int main(int argc, char **argv){

  long b = 2000;
  for (long dx = 5; dx < 40; dx += 5)
    for (long dy = 5; dy < 40; dy += 5){
      do_test(b, dx, dy, stoi(argv[1]), stoi(argv[2]));
      exit(0);
    }
  
}
