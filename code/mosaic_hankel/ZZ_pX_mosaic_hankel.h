#ifndef __ZZ_PX_MOSAIC_HANKEL_H
#define __ZZ_PX_MOSAIC_HANKEL_H

#include <NTL/ZZ_pX.h>
#include <NTL/matrix.h>

#include "lzz_pX_CRT.h"
#include "lzz_p_cauchy_geometric.h" 

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Hankel matrices                                    */
/* stored as                                          */
/*       a5 a4 a3 a2                                  */
/*       a4 a3 a2 a1                                  */
/*       a3 a2 a1 a0                                  */
/*----------------------------------------------------*/
/*----------------------------------------------------*/
class ZZ_p_hankel{

  // n rows, m columns
  long n, m;
  Vec<ZZ_p> data;

 public:
  Vec<ZZ_p> data_rev;

  ZZ_p_hankel(){
    data.SetLength(0);
    n = m = 0;
  }

  /*----------------------------------------------------*/
  /* input vector is as showed above                    */
  /*----------------------------------------------------*/
  ZZ_p_hankel(Vec<ZZ_p>& input, long rows, long cols){
    n = rows;
    m = cols;
    data = input;
    data_rev.SetLength(n+m-1);
    ZZ_pX data_X;
    data_X.rep.SetLength(n+m-1);

    for (long i = 0;i < n+m-1; i++){
      data_rev[i] = input[n+m-2-i];
      data_X.rep[i] = data_rev[i];
    }
    data_X.normalize();
  }

  /*----------------------------------------------------*/
  /* dimensions                                         */
  /*----------------------------------------------------*/
  long NumCols() const{
    return m;
  }
  long NumRows() const{
    return n;
  }
  const ZZ_p& operator ()(long i, long j) const{
    return data[m+n-2-i-j];
  }

};

/*----------------------------------------------------*/
/* turns M into a dense matrix                        */
/*----------------------------------------------------*/
void to_dense(Mat<ZZ_p>& Mdense, const ZZ_p_hankel& M);

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Mosaic Hankel matrices:                            */
/* block matrix where each block is Hankel            */
/*----------------------------------------------------*/
/*----------------------------------------------------*/
class ZZ_p_mosaic_hankel{
  long n, m; // numbers of rows / colums
  long nb, mb; // number of blocks in rows / columns

 public:

  // data[i][j] is the block at i-th row, j-th column
  Vec<Vec<ZZ_p_hankel>> data;

  /*----------------------------------------------------*/
  /* dummy constructor                                  */
  /*----------------------------------------------------*/
  ZZ_p_mosaic_hankel(){
    data.SetLength(0);
    n = m = nb = mb = 0;
  }

  /*----------------------------------------------------*/
  /* copies all data                                    */
  /*----------------------------------------------------*/
  ZZ_p_mosaic_hankel(Vec<Vec<ZZ_p_hankel>> init){
    data = init;
    nb = init.length();
    mb = init[0].length();

    n = 0;
    m = 0;

    for(long i = 0; i < nb; i++)
      n += init[i][0].NumRows();
    for (long j = 0; j < mb; j++)
      m += init[0][j].NumCols();
  }

  /*----------------------------------------------------*/
  /* getters                                            */
  /*----------------------------------------------------*/
  long NumRows() const{
    return n;
  }
  long NumCols() const{
    return m;
  }
  long NumBlockRows() const{
    return nb;
  }
  long NumBlockCols() const{
    return mb;
  }
  long NumRows_of_block(long i) const{
    return data[i][0].NumRows();
  }
  long NumCols_of_block(long i) const{
    return data[0][i].NumCols();
  }

  const ZZ_p& operator ()(long i, long j) const{

    if (i >= n || j >= m)
      Error("matrix indices out of bounds\n");
    long idx, jdx;

    idx = 0;
    while (i >= data[idx][0].NumRows()){
      i -= data[idx][0].NumRows();
      idx++;
    }

    jdx = 0;
    while (j >= data[0][jdx].NumCols()){
      j -= data[0][jdx].NumCols();
      jdx++;
    }
    return data[idx][jdx](i, j);
  }

};

/*----------------------------------------------------*/
/* turns M into a dense matrix                        */
/*----------------------------------------------------*/
void to_dense(Mat<ZZ_p>& Mdense, const ZZ_p_mosaic_hankel& M);

/*----------------------------------------------------*/
/* access to particular rows and columns              */
/*----------------------------------------------------*/
void first_column_of_block(Vec<ZZ_p>& res, long i, const ZZ_p_mosaic_hankel& M);
void last_column_of_block(Vec<ZZ_p>& res, long i, const ZZ_p_mosaic_hankel& M);
void first_row_of_block(Vec<ZZ_p>& res, long i, const ZZ_p_mosaic_hankel& M);
void last_row_of_block(Vec<ZZ_p>& res, long i, const ZZ_p_mosaic_hankel& M);

/*----------------------------------------------------*/
/* G, H such that Z1 M - Z0^t M = G H^t               */
/*----------------------------------------------------*/
void generators(Mat<ZZ_p>& G, Mat<ZZ_p>& H, const ZZ_p_mosaic_hankel& M);

#endif
