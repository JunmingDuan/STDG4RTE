#include "DGFEMSpace1D.h"

void DGFEMSpace1D::solve_leqn(const EMAT& A, const EVEC& rhs, EVEC& u) {
  solver.compute(A);
  u = solver.solve(rhs);
}

void DGFEMSpace1D::STDG(const SOL& I, const double t, const double dt,
      const u_int cell, const u_int m,
    func BDfunc, func q, ST_ele& W) {
  double hi = mesh[cell+1] - mesh[cell];
  if(mu[m] > 0) {
    A = (hi/c)*Kt + (dt*mu[m])*Kxl + (sigma_t*hi*dt)*MM;//all
    rhs = (hi/c)*Ft*I[cell][m] + source(q, cell, m, t, dt);
    rhs += (mu[m]*dt/2.)*STDG_BD(cell, m, t, dt, BDfunc);
  }
  else {
    A = (hi/c)*Kt + (dt*mu[m])*Kxr + (sigma_t*hi*dt)*MM;//all
    rhs = (hi/c)*Ft*I[cell][m] + source(q, cell, m, t, dt);
    rhs -= (mu[m]*dt/2.)*STDG_BD(cell, m, t, dt, BDfunc);
  }
  solve_leqn(A, rhs, W);
}

void DGFEMSpace1D::STDG(const SOL& I, const double t, const double dt,
    const u_int cell, const u_int m,
    const ST_ele& W_uw, func q, ST_ele& W) {
  double hi = mesh[cell+1] - mesh[cell];
  if(mu[m] > 0) {
    A = (hi/c)*Kt + (dt*mu[m])*Kxl + (sigma_t*hi*dt)*MM;//all
    rhs = (hi/c)*Ft*I[cell][m] + (dt*mu[m])*Fxl*W_uw;
    rhs += source(q, cell, m, t, dt);
  }
  else {
    A = (hi/c)*Kt + (dt*mu[m])*Kxr + (sigma_t*hi*dt)*MM;//all
    rhs = (hi/c)*Ft*I[cell][m] + (dt*mu[m])*Fxr*W_uw;
    rhs += source(q, cell, m, t, dt);
  }
  solve_leqn(A, rhs, W);
}


EVEC DGFEMSpace1D::source(func q, u_int cell, u_int m,
    const double t, const double dt) {
  EVEC s(K*K);
  s.setZero();
  double hi = mesh[cell+1] - mesh[cell];
  u_int px, pt;
  VEC<double> x0 = TemQuad_Gauss.points();
  VEC<double> wg = QUADINFO[cell].weight();
  VEC<double> xg = QUADINFO[cell].points();
  VEC<double> lv(2), gv(2);
  lv[0] = -1, lv[1] = 1;
  gv[0] = t, gv[1] = t+dt;

  double tg;
  for(u_int p = 0; p < K*K; ++p) {
    st2coe(p, px, pt);
    local_to_global(x0[pt], lv, gv, &tg);
    s[p] += wg[px]*wg[pt]*q(mu[m], xg[px], tg);
  }
  return s*(hi*dt/4.0);
}

void DGFEMSpace1D::update_sol(const VEC<ST_ele>& W, SOL& I) {
  u_int px, pt;
  for(u_int m = 0; m < M; ++m) {//for each light direction
    for(u_int i = 0; i < Nx; ++i) {//for each cell, update I_{m,i}
      I[i][m].setZero();
      for(u_int p = 0; p < K*K; ++p) {
        st2coe(p, px, pt);
        I[i][m][px] += W_list[i][p]*Lagrange_Poly(1)[pt];
      }
    }//end i
  }
}

double DGFEMSpace1D::RAD_BE_unsteady(const SOL& In, SOL& I, func q,
    const double t, const double dt, func BL, func BR) {//period
  for(u_int m = 0; m < M; ++m) {//for each light direction
    if(mu[m] >= 0) {

      for(u_int i = 0; i < Nx; ++i) {//predictor
        if(i == 0) {
          STDG(I, t, dt, i, m, BL, q, W_list[i]);
        }
        else {
          STDG(I, t, dt, i, m, W_list[i-1], q, W_list[i]);
        }
      }

    }//end ->
    else {
      u_int i = Nx-1;
      while(i >= 0) {
        if(i == Nx-1) {
          STDG(I, t, dt, i, m, BR, q, W_list[i]);
        }
        else {
          STDG(I, t, dt, i, m, W_list[i+1], q, W_list[i]);
        }

        if(i-- == 0) break;
      }//end i
    }//end <-
  }//end M

  update_sol(W_list, I);

  return 1;
}

void DGFEMSpace1D::run_unsteady(func q, func BL, func BR, double t_end) {
  std::cout.precision(16);
  std::cout << std::showpos;
  std::cout.setf(std::ios::scientific);
  double t(0), dt(0);
  //print trouble cells
  std::stringstream s;
  s << "ex1_Nx" << Nx << "_K" << K << "_TC.dat";
  std::string filename(s.str());
  std::ofstream out(filename.c_str());

  while (t < t_end) {
    dt = cal_dt(I, t);
    dt = std::min(dt, t_end-t);
    RAD_BE_unsteady(In, I, q, t, dt, BL, BR);
    Reconstruct_TC(I, In, out);
    I = In;
    //In = I;
    t += dt;
    std::cout << "t: " << t << ", dt: " << dt << std::endl;
  }

  out.close();
}

