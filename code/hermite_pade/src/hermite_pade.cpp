
#include "ZZ_pX_extra.h"
#include "magma_output.h"
#include "ZZ_hermite_pade.h"
#include "vec_ZZ_p_extra.h"
#include "lzz_pXY.h"
#include "lzz_p_extra.h"

NTL_CLIENT

/*----------------------------------------------------------------*/
/* applies a block reversal to v                                  */
/* e.g., type = [1,2] v = [v1,v2,v3,v4,v5]                        */
/* returns [v2,v1,v4,v3,v2] (blocks have length type[i]+1)        */                                    
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::flip_on_type (const Vec<ZZ_p> &v){
  Vec<ZZ_p> r;
  r.SetMaxLength(v.length());
  long acc = 0;
  for (long i = 0; i < type.length(); i++){
    for(long j = type[i]; j >= 0; j--)
      r.append(v[j+acc]);
    acc += type[i] + 1;
  }
  return r;
}

/*----------------------------------------------------------------*/
/* applies a block reversal to v                                  */
/* e.g., type = [1,2] v = [v1,v2,v3,v4,v5]                        */
/* returns [v2,v1,v4,v3,v2] (blocks have length type[i]+1)        */                                    
/*----------------------------------------------------------------*/
Vec<Vec<ZZ>> hermite_pade::flip_on_type (const Vec<Vec<ZZ>> &v){
  Vec<Vec<ZZ>> r;
  r.SetMaxLength(v.length());
  long acc = 0;
  for (long i = 0; i < type.length(); i++){
    for(long j = type[i]; j >= 0; j--)
      r.append(v[j+acc]);
    acc += type[i] + 1;
  }
  return r;
}

/*----------------------------------------------------------------*/
/* splits v according to the current type                         */
/*----------------------------------------------------------------*/
Vec<Vec<ZZ_p>> hermite_pade::split_on_type(const Vec<ZZ_p> &v){
  long acc = 0;
  Vec<Vec<ZZ_p>> result;
  result.SetLength(type.length());
  for (long i = 0; i < type.length(); i++){
    Vec<ZZ_p> f;
    f.SetLength(type[i]+1);
    for (long j = 0; j < type[i]+1; j++)
      f[j] = v[acc+j];
    acc += type[i] + 1;
    result[i] = f;
  }
  return result;
}

/*----------------------------------------------------------------*/
/* splits v according to the current type                         */
/*----------------------------------------------------------------*/
Vec<zz_pX> hermite_pade::split_on_type(const Vec<zz_p> &v){
  long acc = 0;
  Vec<zz_pX> result;
  for (long i = 0; i < type.length(); i++){
    Vec<zz_p> f;
    for (long j = 0; j < type[i]+1; j++)
      f.append(v[acc+j]);
    acc += type[i] + 1;
    result.append(conv<zz_pX>(f));
  }
  return result;
}

/*----------------------------------------------------------------*/
/* splits v according to the current type                         */
/*----------------------------------------------------------------*/
Vec<ZZX> hermite_pade::split_on_type(const Vec<ZZ> &v){
  long acc = 0;
  Vec<ZZX> result;
  for (long i = 0; i < type.length(); i++){
    Vec<ZZ> f;
    for (long j = 0; j < type[i]+1; j++)
      f.append(v[acc+j]);
    acc += type[i] + 1;
    result.append(conv<ZZX>(f));
  }
  return result;
}

/*----------------------------------------------------------------*/
/* switches ZZ_p to ZZ mod p^(2^n)                                */
/* does not save the current context                              */
/*----------------------------------------------------------------*/
void hermite_pade::switch_context(long n){
  if (n < vec_M.length()){ // it has already been computed
    ZZ_p::init(p_powers[n]);
    level = n;
  } 
  else{
    if (n != vec_M.length())
      throw "length mismatch for table of bmc";

    // calculating the new power of p
    ZZ p_new(p);
    long pow2 = 1L << n;
    power(p_new, p_new, pow2); // 2^n isn't going to be very large
    p_powers.append(p_new);
    ZZ_p::init(p_new);
    
    // creating the new bivariate modular comp
    set_up_bmc();
    
    // computing w mod p^(2^n)
    ZZ new_w;
    lift_root_of_unity(new_w, this->w, order, p, pow2);
    ZZ_p w_p, c_p, d_p;
    conv(w_p, new_w);
    conv(c_p, c);
    conv(d_p, d);
    ZZ_pX_Multipoint_FFT X_new (w_p, c_p, sizeX);
    ZZ_pX_Multipoint_FFT Y_new (w_p, d_p, sizeY);
    // update
    
    vec_X_int.append(X_new);
    vec_Y_int.append(Y_new);
    level = n;
  }
}

/*----------------------------------------------------------------*/
/* calls create_bmc and appends the result to vec_M               */
/*----------------------------------------------------------------*/
void hermite_pade::set_up_bmc(){
  vec_M.append(create_bmc());
}

/*----------------------------------------------------------------*/
/* multiplies b by the matrix Y_int^t D_f                         */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mul_Y_right(const Vec<ZZ_p> &b){
  Vec<ZZ_p> x;
  Vec<ZZ_p> f = conv<Vec<ZZ_p>>(this->f); // TODO: e and f should be called e_ZZ and f_ZZ?
  mul_diagonal(x, f, b);
  vec_Y_int[level].mul_left(x, x);
  return x;
}

/*----------------------------------------------------------------*/
/* multiplies b by the matrix D_e X_int                           */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mul_X_right(Vec<ZZ_p> b){
  vec_X_int[level].mul_right(b, b);
  Vec<ZZ_p> e = conv<Vec<ZZ_p>>(this->e);
  Vec<ZZ_p> x;
  mul_diagonal(x, e, b);
  return x;
}

/*----------------------------------------------------------------*/
/* multiplies b by the matrix M                                   */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mul_M_right(const Vec<ZZ_p> &b){
  Vec<ZZ_p> flipped = flip_on_type(b);
  return vec_M[level]->mul_right(flipped);
}

/*----------------------------------------------------------------*/
/* multiplies b by the matrix CL =  D_e X_int M Y_int^t D_f       */
/* (CL is cauchy-geometric-like)                                  */
/* b need not have size CL.NumCols()                              */
/*----------------------------------------------------------------*/
Vec<ZZ_p> hermite_pade::mulA_right(const Vec<ZZ_p>& b){
  Vec<ZZ_p> b_loc = b;
  b_loc.SetLength(sizeY, ZZ_p(0)); // pad it
  Vec<ZZ_p> x = mul_Y_right(b_loc);
  x = mul_M_right(x);
  x= mul_X_right(x);
  return x;
}

/*----------------------------------------------------------------*/
/* if Mx = b mod p^(2^{n-1}), updates x so that Mx = b mod p^(2^n)*/
/*----------------------------------------------------------------*/
void hermite_pade::update_solution(Vec<ZZ>& x, const Vec<ZZ_p> &b, long n){
  Vec<ZZ_p> r = mulA_right(conv<Vec<ZZ_p>>(x)) - b; // mod p^(2^n)
  Vec<ZZ> r_ZZ = conv<Vec<ZZ>>(r);
  for (long i = 0; i < r.length(); i++)
    r_ZZ[i] = r_ZZ[i] / p_powers[n-1];
  Vec<ZZ> x_1;
  DAC(x_1, r_ZZ, n-1);
  x = x - p_powers[n-1] * x_1;
}

/*----------------------------------------------------------------*/
/* solves for Mx = b mod p^(2^n)                                  */
/*----------------------------------------------------------------*/
void hermite_pade::DAC(Vec<ZZ> &x, const Vec<ZZ>& b_in, long n){

  if (n == 0){ // we are mod p
    zz_pContext push;
    ctx.restore();
    Vec<zz_p> x_zz_p;
    Vec<zz_p> b_zz_p = conv<Vec<zz_p>>(b_in);
    b_zz_p.SetLength(rank); // TODO: why do we need this?
    invA.mul_right(x_zz_p, b_zz_p);
    x = conv<Vec<ZZ>>(x_zz_p);
    return;
  }

  long old_n = level;
  switch_context(n); 

  Vec<ZZ_p> b_n = conv<Vec<ZZ_p>>(b_in);
  Vec<ZZ> b = conv<Vec<ZZ>>(b_n); // this reduces b_in mod p^(2^n).
  DAC(x, b, n-1); 

  update_solution(x, b_n, n);
  switch_context(old_n); 
}

/*----------------------------------------------------------------*/
/* checks if every entry can be reconstructed                     */
/* if so, check whether the solution cancels the system mod p2    */
/*----------------------------------------------------------------*/
bool hermite_pade::reconstruct_and_check(Vec<ZZX> & sol_poly, const Vec<ZZ_p> &v, long n){
  Vec<Vec<ZZ>> sol;
  sol.SetLength(0);

  
  ZZ bound;
  if (n == 0)
    bound = SqrRoot(ZZ(p));
  else
    bound = p_powers[n-1];

  ZZ_p denom = to_ZZ_p(0);
  for (long i = 0; i < v.length(); i++)
    if (v[i] != 0){
      denom = v[i];
      break;
    }
  ZZ_p i_den = 1 / denom;

  for (long i = 0; i < v.length(); i++){
    ZZ a,b;
    try{
      long result = ReconstructRational(a, b, conv<ZZ>(v[i] * i_den), p_powers[n], bound, bound);
      if (result == 0) 
	return false;
      Vec<ZZ> temp;
      temp.append(a);
      temp.append(b);
      sol.append(temp);
    }
    catch(...){
      return false;
    }
  }

  sol = flip_on_type(sol);

  Vec<ZZ> sol_ZZ;
  sol_ZZ.SetLength(sol.length());

  ZZ ell{1};
  for (long i = 0; i < sol.length(); i++)
    ell = (sol[i][1] * ell) / GCD(sol[i][1], ell);
  for (long i = 0; i < sol.length(); i++)
    sol_ZZ[i] = (sol[i][0] * ell) / sol[i][1];

#if true  
  static Vec<ZZ> sol_ZZ_old;
  if (sol_ZZ != sol_ZZ_old){
    sol_ZZ_old = sol_ZZ;
    return false;
  }
#else
  ZZ_pContext push;
  ctx2.restore();
  
  Vec<ZZ_p> sol_ZZ_p;
  for (long int i = 0; i < sol.length(); i++)
    sol_ZZ_p.append(conv<ZZ_p>(sol_ZZ[i]));

  auto bmc = create_bmc();
  auto x = bmc->mul_right(sol_ZZ_p);
  for (long int i = 0; i < x.length(); i++)
    if (x[i] != ZZ_p(0))
      return false;
#endif

  sol_poly = split_on_type(sol_ZZ);

  return true;
}

/*----------------------------------------------------------------*/
/* computes a random solution of the system modulo p              */
/*----------------------------------------------------------------*/
void hermite_pade::random_solution_mod_p(Vec<zz_pX> & pol){
  pol.SetLength(type.length());

  Vec<ZZ_p> ex1; // mult with A to get a column
  Vec<ZZ_p> add1;
  Vec<ZZ> sol1;

  ex1.SetLength(sizeY);
  add1.SetLength(sizeY - rank);
  for (long i = 0; i < rank; i++)
    ex1[i] = 0;
  for (long i = 0; i < add1.length(); i++){
    add1[i] = random_ZZ_p();
    ex1[rank+i] = add1[i];
  }
  ex1 = mulA_right(ex1);

  DAC(sol1, conv<Vec<ZZ>>(ex1), 0);
  conv(ex1, sol1);
  ex1.SetLength(sizeY, ZZ_p(0));
  for (long i = 0; i < add1.length(); i++)
    ex1[rank+i] = -add1[i];

  ex1 = flip_on_type(mul_Y_right(ex1));
  Vec<zz_p> v = conv<Vec<zz_p>>(conv<Vec<ZZ>>(ex1));

  zz_p denom = to_zz_p(0);
  for (long i = 0; i < v.length(); i++)
    if (v[i] != 0){
      denom = v[i];
      break;
    }
  if (denom == 0)
    return;
  for (long i = 0; i < v.length(); i++)
    v[i] /= denom;
  pol = split_on_type(v);
}

/*----------------------------------------------------------------*/
/* computes a solution of the system modulo p                     */
/* assumes that the nullity is one                                */
/*----------------------------------------------------------------*/
void hermite_pade::random_solution(Vec<ZZX> &sol_poly){
  if ((NumCols() - Rank()) != 1)
    throw "lifting works only in case of nullity 1";

  zz_pContext zz_p_push;
  zz_pContext ZZ_p_push;

  ctx.restore();
  long n = 0; // start at p^2^n
  switch_context(n);


  Vec<ZZ_p> extractor; // mult with A to get a column
  extractor.SetLength(sizeY, ZZ_p(0));
  extractor[rank] = 1; // just for now, take the last column

  Vec<ZZ_p> b = mulA_right(extractor); // b is the last column of A
  Vec<ZZ> x_ZZ, b_ZZ;
  b_ZZ = conv<Vec<ZZ>> (b);
  DAC(x_ZZ, b_ZZ, n); // solution mod p

  Vec<ZZ_p> x = conv<Vec<ZZ_p>> (x_ZZ);
  x.SetLength(sizeY, ZZ_p(0));
  x[rank] = -1;
  Vec<ZZ_p> soln = mul_Y_right(x);


  // loop until we get enough prec
  while (!reconstruct_and_check(sol_poly, soln, n)){
    switch_context(++n);

    Vec<ZZ_p> extractor; // mult with A to get a column
    extractor.SetLength(sizeY, ZZ_p(0));
    extractor[rank] = 1; // just for now, take the last column
    Vec<ZZ_p> b = mulA_right(extractor); // b is the last column of A

    update_solution(x_ZZ, b, n);

    Vec<ZZ_p> x = conv<Vec<ZZ_p>> (x_ZZ);
    x.SetLength(sizeY,ZZ_p(0));
    x[rank] = -1;
    soln.kill();
    soln = mul_Y_right(x);

    ZZ_p first;
    for (long i = 0; i < soln.length(); i++)
      if (soln[i]._ZZ_p__rep % p_powers[0] != ZZ(0)){
        first = 1/soln[i];
        break;
      }
    soln *= first;
  }
}


/*----------------------------------------------------------------*/
/* Rank                                                           */
/*----------------------------------------------------------------*/
long hermite_pade::Rank() const {
  return rank;
}

/*----------------------------------------------------------------*/
/* Dimensions: number of rows                                     */
/*----------------------------------------------------------------*/
long hermite_pade::NumRows() const {
  return sizeX;
}

/*----------------------------------------------------------------*/
/* Dimensions: number of columns                                  */
/*----------------------------------------------------------------*/
long hermite_pade::NumCols() const {
  return sizeY;
}

/*----------------------------------------------------------------*/
/* sets up the field, contexts, ..                                */
/*----------------------------------------------------------------*/
hermite_pade::hermite_pade(long fft_index) {
  level = 0;
  ctx2 = ZZ_pContext(ZZ(9001));
  zz_p::FFTInit(fft_index);
  ctx = zz_pContext(INIT_FFT, fft_index);
  p = zz_p::modulus();
  ZZ_p::init(ZZ(p));
  p_powers.append(ZZ(p));
}
