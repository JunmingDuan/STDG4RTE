#include "DGFEMSpace1D.h"

void DGFEMSpace1D::Boundary(func BL, func BR, double& BD_L, double& BD_R, const double mu, const double t, const double dt) {
  VEC<VEC<double> > pw = QUADINFO[0].LG(K);
  VEC<double> gv(2), lv(2);
  lv[0] = -1, lv[1] = 1;
  gv[0] = t, gv[1] = dt+t;
  VEC<double> gt(K);
  for(u_int j = 0; j < K; ++j) {
    local_to_global(pw[0][j], lv, gv, &gt[j]);
  }

  BD_L = 0; BD_R = 0;
  for(u_int g = 0; g < K; ++g) {
    BD_L += pw[1][g]*BL(mu, mesh[0], gt[g]);
    BD_R += pw[1][g]*BR(mu, mesh[Nx], gt[g]);
  }
  BD_L *= mu*dt/2.;
  BD_R *= mu*dt/2.;

}

EVEC DGFEMSpace1D::STDG_BD(const u_int cell, const u_int m, const double t, const double dt, func BDfunc) {
  VEC<VEC<double> > pw = QUADINFO[0].LG(K);
  VEC<double> gv(2), lv(2);
  lv[0] = -1, lv[1] = 1;
  gv[0] = t, gv[1] = dt+t;
  EVEC flux(K*K);
  flux.setZero();
  u_int px, pt;
  for(u_int p = 0; p < K*K; ++p) {
    double tg;
    st2coe(p, px, pt);
    local_to_global(pw[0][pt], lv, gv, &tg);
    if(cell == 0)
      flux[p] = pw[1][pt]*BDfunc(mu[m], mesh[0], tg)*Lagrange_Poly(-1)[px];
    else if(cell == Nx-1)
      flux[p] = pw[1][pt]*BDfunc(mu[m], mesh[Nx], tg)*Lagrange_Poly(1)[px];
  }
  return flux;
}

