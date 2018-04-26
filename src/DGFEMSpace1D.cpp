/**
 * @file DGFEMSpace1D.cpp
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0 for 1D RTE
 * @
 * @date 2018-04-26
 */

#include <stdio.h>
#include <stdlib.h>
#include "DGFEMSpace1D.h"

DGFEMSpace1D::DGFEMSpace1D(u_int Nx, double xl, double xr)
  : Nx(Nx), xl(xl), xr(xr) {
  mesh.resize(Nx+1);
  h = (xr - xl)/Nx;
  for(u_int i = 0; i < Nx+1; ++i) {
    mesh[i] = i*h;
  }
  mu.resize(M);
  wgt.resize(M);
  I.resize(Nx); I1.resize(Nx);
  In.resize(Nx); I_last.resize(Nx);
  I_av.resize(Nx);
  W_list.resize(Nx);
  W_list_last.resize(Nx);
  for(u_int i = 0; i < Nx; ++i) {
    I[i].resize(M); I1[i].resize(M);
    In[i].resize(M); I_last[i].resize(M);
    I_av[i].resize(M);
    W_list[i].resize(K*K);
    W_list_last[i].resize(K*K);
    for(u_int m = 0; m < M; ++m) {
      I[i][m].resize(K); I1[i][m].resize(K);
      In[i][m].resize(K); I_last[i][m].resize(K);
    }
  }
  Kt.resize(K*K, K*K);
  Kxl.resize(K*K, K*K);
  Fxl.resize(K*K, K*K);
  Kxr.resize(K*K, K*K);
  Fxr.resize(K*K, K*K);
  MM.resize(K*K, K*K);
  Ft.resize(K*K, K);
  A.resize(K, K);
  rhs.resize(K);
}

void DGFEMSpace1D::BuildQuad(u_int np) {
  VEC<VEC<double> > pw;
  //Legendre Gauss points
  pw = QUADINFO[0].LG(np);
  TemQuad_Gauss.set_np(np);
  TemQuad_Gauss.set_points(pw[0]);
  TemQuad_Gauss.set_weight(pw[1]);
  QUADINFO.resize(Nx);
  VEC<double> gv(2), lv(2);
  lv[0] = -1, lv[1] = 1;
  for(u_int i = 0; i < Nx; ++i) {
    QUADINFO[i].set_np(np);
    gv[0] = mesh[i], gv[1] = mesh[i+1];
    QUADINFO[i].set_weight(pw[1]);
    VEC<double> x(np);
    for(u_int j = 0; j < np; ++j) {
      local_to_global(pw[0][j], lv, gv, &x[j]);
    }
    QUADINFO[i].set_points(x);
    QUADINFO[i].set_jacobi( local_to_global_jacobian(lv, gv),
        global_to_local_jacobian(lv, gv) );
  }
  //direction
  sum_wm = 2;
  if(M > 1) {
    pw = TemQuad_Gauss.LG(M);
    for(u_int m = 0; m < M; ++m) {
      mu[m] = pw[0][m];
      wgt[m] = pw[1][m];
    }
  }
  else {
    mu[0] = 1;
    wgt[0] = 2;
  }

}

void DGFEMSpace1D::st2coe(const u_int q, u_int& qx, u_int& qt) {
  qx = q%K;
  qt = (q-qx)/K;
}

void DGFEMSpace1D::BuildUniMAT() {
  u_int qx, qt, px, pt;
  VEC<double> x = TemQuad_Gauss.points();
  VEC<double> w = TemQuad_Gauss.weight();
  EVEC pl = Lagrange_Poly(-1);
  EVEC pr = Lagrange_Poly(1);
  EVEC prime;
  Kt.setZero();
  Kxl.setZero();
  Fxl.setZero();
  Kxr.setZero();
  Fxr.setZero();
  MM.setZero();
  Ft.setZero();
  for(u_int p = 0; p < K*K; ++p) {
    for(u_int q = 0; q < K*K; ++q) {
      st2coe(p, px, pt);
      st2coe(q, qx, qt);
      if(pt == qt) {
        prime = Lagrange_PolyG(x[qx]);
        Kxl(p,q) = 0.5*w[pt]*pr[px]*pr[qx] -
          0.5*w[qx]*w[qt]*prime[px];
        Fxl(p,q) = 0.5*w[pt]*pl[px]*pr[qx];
        Kxr(p,q) = - 0.5*w[pt]*pl[px]*pl[qx] -
          0.5*w[qx]*w[qt]*prime[px];
        Fxr(p,q) = - 0.5*w[pt]*pr[px]*pl[qx];
      }
      if(px == qx) {
        prime = Lagrange_PolyG(x[qt]);
        Kt(p,q) = 0.5*w[px]*pr[pt]*pr[qt]
          - 0.5*w[qx]*w[qt]*prime[pt];
      }
      if(p == q) {
        MM(p,q) = 0.25*w[px]*w[pt];
      }
    }
  }
  for(u_int p = 0; p < K*K; ++p) {
    for(u_int l = 0; l < K; ++l) {
      st2coe(p, px, pt);
      if(l == px) {
        Ft(p,l) = 0.5*w[l]*pl[pt];
      }
    }
  }

  //std::cout << "Fx: " << Fx << std::endl;
  //abort();
}

void DGFEMSpace1D::Projection(u_int cell, func I0, double t, bU& u) {
  VEC<double> p = QUADINFO[cell].points();
  for(u_int m = 0; m < M; ++m) {
    for(u_int g = 0; g < K; ++g) {
      u[m][g] = I0(mu[m], p[g], t);
    }
  }
}

void DGFEMSpace1D::DG2av(const SOL& I, VEC<VEC<double>>& I_av) {
  VEC<double> w = TemQuad_Gauss.weight();
  for(u_int m = 0; m < M; ++m) {
    for(u_int i = 0; i < Nx; ++i) {
      I_av[i][m] = 0;
      for(u_int g = 0; g < K; ++g) {
        I_av[i][m] += w[g]*I[i][m][g];
      }
      I_av[i][m] /= 2.0;
    }
  }
}

EVEC DGFEMSpace1D::Composition(const SOL& I, u_int cell, double x) {
  EVEC u(M);
  u.setZero();
  VEC<double> lv(2), gv(2);
  EVEC V;
  lv[0] = -1, lv[1] = 1;
  gv[0] = mesh[cell], gv[1] = mesh[cell+1];

  double lp;
  global_to_local(x, lv, gv, &lp);
  V = Lagrange_Poly(lp);
  for(u_int m = 0; m < M; ++m) {
    for(u_int k = 0; k < K; ++k) {
      u[m] += I[cell][m][k]*V[k];
    }
  }

  return u;
}

void DGFEMSpace1D::init(func I0) {
  for(u_int i = 0; i < Nx; ++i) {
    Projection(i, I0, 0, I[i]);
  }
  In = I;
  I_last = I;
}

double DGFEMSpace1D::cal_dt(const SOL& I, const double t) {
  //return 1e-3;
  double CFL(1);
  return CFL*h;
}

void DGFEMSpace1D::print_solution_integral(const SOL& I, std::ostream& os) {
  os.precision(16);
  os << std::showpos;
  os.setf(std::ios::scientific);
  VEC<double> x = TemQuad_Gauss.points();
  for(u_int i = 0; i < Nx; ++i) {
    VEC<double> w = QUADINFO[i].weight();
    VEC<double> p = QUADINFO[i].points();
    for(u_int g = 0; g < x.size(); ++g) {
      os << p[g] << " "  << w[g];
      for(u_int m = 0; m < M; ++m) {
        os  << " " << I[i][m][g];
      }
      os << "\n";
    }
    os << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
}

void DGFEMSpace1D::print_solution_average(std::ostream& os) {
  os.precision(16);
  os << std::showpos;
  os.setf(std::ios::scientific);
  for(u_int i = 0; i < Nx; ++i) {
    os << 0.5*(mesh[i]+mesh[i+1]) << " " << I[i][0] << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
}

void DGFEMSpace1D::print_sol_ex_integral(std::ostream& os, func I0, const double t) {
  os.precision(16);
  os << std::showpos;
  os.setf(std::ios::scientific);
  I1 = I;
  for(u_int i = 0; i < Nx; ++i) {
    VEC<double> w = QUADINFO[i].weight();
    VEC<double> p = QUADINFO[i].points();
    for(u_int g = 0; g < w.size(); ++g) {
      os << p[g] << " "  << w[g];
      for(u_int m = 0; m < M; ++m) {
        os  << " " << I1[i][m][g];
      }
      for(u_int m = 0; m < M; ++m) {
        os << " " << I0(mu[m], p[g], t);
      }
      os << "\n";
    }
    os << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
}


