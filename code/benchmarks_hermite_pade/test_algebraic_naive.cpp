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
void random_F_and_solution_mod_p(ZZXY & F, ZZX & g, long dx, long dy, long b){
  zz_p::FFTInit(0);
  ZZ p{zz_p::modulus()};
  while (p < (ZZ{1} << (2*b)))
    p = p*p;
  p = p*p;
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
/* single test. mode_lifting is as in ZZ_hermite_pade.h       */
/*------------------------------------------------------------*/
void do_test(long b, long dx, long dy, long mode_lifting){

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
  ZZ p_back = ZZ_p::modulus();
  
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

  ZZ_p::init(p_back);
  Vec<ZZX> vecF;
  ZZX tmp;
  SetCoeff(tmp, 0, ZZ{1});
  for (long i = 0; i < type.length(); i++){
    vecF.append(tmp);
    tmp = conv<ZZX>(conv<ZZ_pX>(trunc(tmp * f, sigma)));
  }

  hermite_pade_general hpg = hermite_pade_general(vecF, type, sigma, witness, 0);
  
  hpg.switch_mode(mode_lifting);

  cout << "dx=" << dx << " dy=" << dy << endl;
  cout << "type " << type << endl;
  cout << "dimensions:= " << hpg.NumRows() << " " << hpg.NumCols() << endl;
  cout << "nullity " << hpg.NumCols() - hpg.Rank() << endl;

  Vec<ZZX> sol;
  hpg.random_solution(sol);
  cout << "first coefficient of sol: " << coeff(sol[0], 0) << endl;
  cout << "first coefficient of F: " << coeff(F.coeffX[0], 0) << endl;
 }

/*------------------------------------------------------------*/
/* usage:                                                     */
/*------------------------------------------------------------*/
int main(int argc, char **argv){

  long b = 20;
  for (long dx = 15; dx < 40; dx += 5)
    for (long dy = 15; dy < 40; dy += 5){
      do_test(b, dx, dy, stoi(argv[1]));
      exit(0);
    }
  
}
