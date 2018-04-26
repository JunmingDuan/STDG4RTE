#include "DGFEMSpace1D.h"

double DGFEMSpace1D::cal_res(const SOL& s1, const SOL& s2) {
  double res(0);
  for(u_int i = 0; i < Nx; ++i) {
    for(u_int d = 0; d < M; ++d) {
      for(u_int k = 0; k < K; ++k) {
        res = std::max(fabs(s1[i][d][k]-s2[i][d][k]), res);
      }
    }
  }
  return res;
}

double DGFEMSpace1D::cal_res(const VEC<ST_ele>& w1, const VEC<ST_ele>& w2) {
  double res(0);
  for(u_int i = 0; i < Nx; ++i) {
    for(u_int k = 0; k < K*K; ++k) {
      res = std::max(fabs(w1[i][k]-w2[i][k]), res);
    }
  }
  return res;
}


